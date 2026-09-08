"""Tests for the striding parallel driver.

The contract is narrow and worth stating plainly: dividing the work must not
change the answer.  Every integer accumulator and the whole coverage vector
must come back identical however many workers ran, because all of them are
exact sums over a partition of the fragments.  ``p_seq_true_sum`` is the one
exception, and PARALLEL_SUM in pytransrate.read_metrics says why.
"""

from __future__ import annotations

import multiprocessing

import numpy as np
import pysam
import pytest

from pytransrate.assign import assign_fragments, group_by_fragment
from pytransrate.bam_metrics import (
    MalformedRecordStats,
    accumulate_metrics,
    estimate_realistic_distance,
    iter_alignments,
)
from pytransrate.read_metrics import _accumulate_parallel, _worker_count

REFS = [("contigA", 600), ("contigB", 500), ("contigC", 400)]

#: Accumulators that must survive splitting exactly, not approximately.
EXACT_FIELDS = (
    "reads_mapped",
    "bases_mapped",
    "fragments_mapped",
    "both_mapped",
    "properpair",
    "bridges",
    "good",
    "clipped_alignments",
    "clipped_bases",
    "leading_clipped_bases",
    "bases_uncovered",
)

pytestmark = pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(),
    reason="the shared accumulators are mmaps inherited across fork",
)


def _pair(name, ref_id, start, mate_start, *, nm=0, cigar="100M", secondary=False):
    """One read-1/read-2 pair, as an aligner emits it: mates adjacent."""
    out = []
    for read1 in (True, False):
        read = pysam.AlignedSegment()
        read.query_name = name
        flag = 1 | 2 | (64 if read1 else 128) | (16 if not read1 else 0)
        flag |= 32 if read1 else 0
        if secondary:
            flag |= 256
        read.flag = flag
        read.reference_id = ref_id
        read.reference_start = start if read1 else mate_start
        read.next_reference_id = ref_id
        read.next_reference_start = mate_start if read1 else start
        read.mapping_quality = 60
        read.cigarstring = cigar
        read.query_sequence = "A" * 100
        read.query_qualities = pysam.qualitystring_to_array("I" * 100)
        read.set_tag("NM", nm, value_type="i")
        out.append(read)
    return out


def _library(tmp_path, n_fragments=60, multimap=True):
    """A read-ordered BAM with multi-mapping fragments, as snap emits."""
    reads = []
    for index in range(n_fragments):
        ref = index % len(REFS)
        start = 10 + (index % 7) * 30
        reads += _pair(f"frag{index}", ref, start, start + 150, nm=index % 5)
        if multimap and index % 3 == 0:
            other = (ref + 1) % len(REFS)
            reads += _pair(
                f"frag{index}", other, start, start + 150, nm=4, secondary=True
            )
    path = str(tmp_path / "reads.bam")
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": name, "LN": length} for name, length in REFS],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for read in reads:
            out.write(read)
    return path


def _expression():
    return {name: {"eff_count": 10.0 + index} for index, (name, _) in enumerate(REFS)}


def _serial(bam_path):
    names = [name for name, _ in REFS]
    lengths = [length for _, length in REFS]
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return accumulate_metrics(
            names,
            lengths,
            assign_fragments(
                iter_alignments(bam, bam_path, MalformedRecordStats()),
                names,
                _expression(),
            ),
            realistic_distance=estimate_realistic_distance(bam_path),
        )


def _parallel(bam_path, workers, stats=None):
    names = [name for name, _ in REFS]
    lengths = [length for _, length in REFS]
    return _accumulate_parallel(
        bam_path,
        names,
        lengths,
        _expression(),
        realistic_distance=estimate_realistic_distance(bam_path),
        nullprior=0.7,
        workers=workers,
        malformed=stats or MalformedRecordStats(),
    )


# ---------------------------------------------------------------------------
# Splitting the work does not change the answer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("workers", [2, 3, 4, 7])
def test_parallel_matches_serial_exactly(tmp_path, workers):
    bam_path = _library(tmp_path)
    want = _serial(bam_path)
    got = _parallel(bam_path, workers)

    assert [c.name for c in got] == [c.name for c in want]
    for expected, actual in zip(want, got):
        for field in EXACT_FIELDS:
            assert getattr(actual, field) == getattr(expected, field), (
                f"{field} on {expected.name} with {workers} workers"
            )
        np.testing.assert_array_equal(
            np.asarray(actual.coverage), np.asarray(expected.coverage)
        )


@pytest.mark.parametrize("workers", [2, 3, 5])
def test_segmentation_is_bit_identical(tmp_path, workers):
    """p_not_segmented comes from integer coverage, so it must not drift."""
    bam_path = _library(tmp_path)
    for expected, actual in zip(_serial(bam_path), _parallel(bam_path, workers)):
        assert actual.p_not_segmented == expected.p_not_segmented


@pytest.mark.parametrize("workers", [2, 3, 5])
def test_p_seq_true_moves_only_in_the_last_bits(tmp_path, workers):
    """The one float sum. See PARALLEL_SUM."""
    bam_path = _library(tmp_path)
    for expected, actual in zip(_serial(bam_path), _parallel(bam_path, workers)):
        assert actual.p_seq_true() == pytest.approx(expected.p_seq_true(), rel=1e-12)


def test_more_workers_than_fragments(tmp_path):
    """Workers that draw no fragments must still contribute a zeroed share."""
    bam_path = _library(tmp_path, n_fragments=3, multimap=False)
    want = _serial(bam_path)
    got = _parallel(bam_path, workers=8)
    for expected, actual in zip(want, got):
        for field in EXACT_FIELDS:
            assert getattr(actual, field) == getattr(expected, field)


def test_a_bam_with_no_alignments(tmp_path):
    bam_path = _library(tmp_path, n_fragments=0)
    for expected, actual in zip(_serial(bam_path), _parallel(bam_path, workers=4)):
        assert actual.reads_mapped == expected.reads_mapped == 0
        assert actual.bases_uncovered == expected.bases_uncovered


# ---------------------------------------------------------------------------
# The stride itself
# ---------------------------------------------------------------------------


def test_strides_partition_the_fragments(tmp_path):
    """Every fragment must be claimed by exactly one worker."""
    bam_path = _library(tmp_path, n_fragments=37)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        whole = [name for name, _ in group_by_fragment(bam.fetch(until_eof=True))]

    for workers in (2, 3, 8):
        seen = []
        for offset in range(workers):
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                seen += [
                    name
                    for name, _ in group_by_fragment(
                        bam.fetch(until_eof=True), stride=workers, offset=offset
                    )
                ]
        assert sorted(seen) == sorted(whole)
        assert len(seen) == len(set(seen))


def test_a_stride_keeps_whole_fragments(tmp_path):
    """A worker must see all of a fragment's records or none of them."""
    bam_path = _library(tmp_path, n_fragments=12)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        whole = {
            name: len(batch)
            for name, batch in group_by_fragment(bam.fetch(until_eof=True))
        }
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for name, batch in group_by_fragment(
            bam.fetch(until_eof=True), stride=3, offset=1
        ):
            assert len(batch) == whole[name]


def test_assignment_is_unaffected_by_the_stride(tmp_path):
    """A fragment's winner depends on the fragment and the priors, not who ran it."""
    bam_path = _library(tmp_path, n_fragments=20)
    names = [name for name, _ in REFS]
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        whole = {
            read.query_name: read.reference_id
            for read in assign_fragments(bam.fetch(until_eof=True), names, _expression())
        }
    strided = {}
    for offset in range(4):
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in assign_fragments(
                bam.fetch(until_eof=True), names, _expression(), stride=4, offset=offset
            ):
                strided[read.query_name] = read.reference_id
    assert strided == whole


# ---------------------------------------------------------------------------
# Failure and configuration
# ---------------------------------------------------------------------------


def test_a_failing_worker_raises_in_the_parent(tmp_path):
    """A worker's exception must not be swallowed into a silent wrong answer."""
    bam_path = _library(tmp_path, n_fragments=5)
    with pytest.raises(FileNotFoundError):
        _parallel(str(tmp_path / "absent.bam"), workers=2)
    # The driver's state must not leak into the next call.
    from pytransrate.read_metrics import _WORKER

    assert _WORKER == {}
    assert _parallel(bam_path, workers=2)


def test_skipped_records_are_counted_once(tmp_path):
    """Every worker reads every record, so only one may report the skips."""
    bam_path = _library(tmp_path, n_fragments=6)
    stats = MalformedRecordStats()
    _parallel(bam_path, workers=4, stats=stats)
    assert stats.skipped == 0


@pytest.mark.parametrize("threads,expected", [(None, 1), (0, 1), (1, 1), (2, 2), (16, 16)])
def test_worker_count(threads, expected):
    assert _worker_count(threads) == expected
