"""Tests for the bam-read replacement.

Coverage is validated against ``samtools depth -a`` -- a genuinely
independent implementation, bundled with pysam -- rather than against
bam-read, which is the artefact being replaced and is wrong about soft clips
(see SOFT_CLIP_FIX in pytransrate.bam_metrics).
"""

from __future__ import annotations

import numpy as np
import pysam
import pytest

from pytransrate.bam_metrics import (
    CSV_COLUMNS,
    ContigMetrics,
    compute_bam_metrics,
    estimate_realistic_distance,
)

REFS = [("contigA", 400), ("contigB", 300)]


def _make_read(
    name,
    ref_id,
    start,
    cigar,
    seq_len,
    *,
    flag=0,
    mate_ref_id=-1,
    mate_start=-1,
    nm=None,
    mapq=60,
):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = flag
    read.reference_id = ref_id
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigarstring = cigar
    read.next_reference_id = mate_ref_id
    read.next_reference_start = mate_start
    read.query_sequence = "A" * seq_len
    read.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
    if nm is not None:
        read.set_tag("NM", nm, value_type="i")
    return read


def _write_bam(path, reads, refs=REFS, sort=True):
    """Write a BAM.

    ``sort=False`` keeps records in the order given, which is what
    ``estimate_realistic_distance`` requires -- it pairs *adjacent* records by
    name, so a coordinate sort would separate the mates.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate" if sort else "unsorted"},
        "SQ": [{"SN": n, "LN": ln} for n, ln in refs],
    }
    if not sort:
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for read in reads:
                out.write(read)
        return str(path)

    unsorted_path = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.sort("-o", str(path), unsorted_path)
    pysam.index(str(path))
    return str(path)


def _samtools_depth(bam_path, refs=REFS):
    """Ground-truth coverage from samtools, as {contig: np.ndarray}.

    ``-J`` is deliberately omitted: it would count bases spanned by a
    deletion as covered, which is neither what pileup.cpp did nor what the
    transrate score means by coverage.
    """
    lengths = dict(refs)
    out = {name: np.zeros(ln, dtype=np.int64) for name, ln in refs}
    text = pysam.depth("-a", "-Q", "0", "-q", "0", bam_path)
    for line in text.splitlines():
        if not line.strip():
            continue
        name, pos, depth = line.split("\t")
        if name in out and int(pos) <= lengths[name]:
            out[name][int(pos) - 1] = int(depth)
    return out


def _coverage_by_name(contigs):
    return {c.name: np.asarray(c.coverage, dtype=np.int64) for c in contigs}


# ---------------------------------------------------------------------------
# Coverage vs. samtools
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cigar,seq_len",
    [
        ("100M", 100),
        ("50M10D50M", 100),
        ("50M10I40M", 100),
        ("10S90M", 100),
        ("90M10S", 100),
        ("10S80M10S", 100),
        ("5S30M5D30M10I20M5S", 100),
        ("40M20N40M", 80),
        ("100=", 100),
        ("50=10X40=", 100),
    ],
)
def test_coverage_matches_samtools(tmp_path, cigar, seq_len):
    bam = _write_bam(
        tmp_path / "one.bam",
        [_make_read("r1", 0, 50, cigar, seq_len, nm=0)],
    )
    contigs = compute_bam_metrics(bam)
    got = _coverage_by_name(contigs)
    want = _samtools_depth(bam)
    np.testing.assert_array_equal(got["contigA"], want["contigA"])
    np.testing.assert_array_equal(got["contigB"], want["contigB"])


def test_coverage_matches_samtools_on_a_pile_of_reads(tmp_path):
    rng = np.random.default_rng(20260906)
    cigars = ["100M", "10S90M", "90M10S", "50M10D50M", "50M10I40M", "20S60M20S"]
    reads = []
    for i in range(300):
        cigar = cigars[i % len(cigars)]
        ref_id = int(rng.integers(0, 2))
        max_start = REFS[ref_id][1] - 120
        reads.append(
            _make_read(
                f"r{i}", ref_id, int(rng.integers(0, max_start)), cigar, 100, nm=0
            )
        )
    bam = _write_bam(tmp_path / "many.bam", reads)

    got = _coverage_by_name(compute_bam_metrics(bam))
    want = _samtools_depth(bam)
    for name, _ in REFS:
        np.testing.assert_array_equal(got[name], want[name])


def test_bases_uncovered_matches_samtools(tmp_path):
    reads = [
        _make_read("r1", 0, 0, "20S100M", 120, nm=0),
        _make_read("r2", 0, 200, "100M20S", 120, nm=0),
    ]
    bam = _write_bam(tmp_path / "gap.bam", reads)
    contigs = compute_bam_metrics(bam)
    want = _samtools_depth(bam)
    by_name = {c.name: c for c in contigs}
    for name, _ in REFS:
        assert by_name[name].bases_uncovered == int(
            np.count_nonzero(want[name] == 0)
        ), name


# ---------------------------------------------------------------------------
# The soft-clip defect specifically
# ---------------------------------------------------------------------------


def _cpp_coverage(read, ref_length):
    """pileup.cpp's addAlignment, including its 'D' || 'S' cursor bug."""
    coverage = np.zeros(ref_length, dtype=np.int64)
    pos = read.reference_start
    for op, op_len in read.cigartuples:
        if op == 0:  # M
            for _ in range(op_len):
                if pos < ref_length:
                    coverage[pos] += 1
                pos += 1
        elif op in (2, 4):  # D or S -- the bug
            pos += op_len
    return coverage


def test_leading_soft_clip_shifted_coverage_in_the_cpp(tmp_path):
    """A leading soft clip displaced coverage by the clip length."""
    read = _make_read("r1", 0, 100, "20S80M", 100, nm=0)
    bam = _write_bam(tmp_path / "clip.bam", [read])

    ours = _coverage_by_name(compute_bam_metrics(bam))["contigA"]
    theirs = _cpp_coverage(read, 400)
    truth = _samtools_depth(bam)["contigA"]

    np.testing.assert_array_equal(ours, truth)
    assert not np.array_equal(theirs, truth)

    # Ours starts at POS; theirs starts 20bp further right.
    assert int(np.flatnonzero(ours)[0]) == 100
    assert int(np.flatnonzero(theirs)[0]) == 120


def test_no_soft_clips_means_no_disagreement(tmp_path):
    """Why this stayed hidden under snap-aligner 1.0dev.96."""
    read = _make_read("r1", 0, 100, "50M10D40M", 90, nm=0)
    bam = _write_bam(tmp_path / "noclip.bam", [read])
    ours = _coverage_by_name(compute_bam_metrics(bam))["contigA"]
    np.testing.assert_array_equal(ours, _cpp_coverage(read, 400))


# ---------------------------------------------------------------------------
# Pair-level counters
# ---------------------------------------------------------------------------

_PAIRED = 1
_PROPER = 2
_MATE_UNMAPPED = 8
_REVERSE = 16
_MATE_REVERSE = 32
_READ1 = 64
_READ2 = 128


def _fr_pair(name, ref_id, pos1, pos2, *, proper=True, mate_ref_id=None):
    """A forward/reverse pair with read 1 to the left."""
    mate_ref = ref_id if mate_ref_id is None else mate_ref_id
    flag1 = _PAIRED | _READ1 | _MATE_REVERSE | (_PROPER if proper else 0)
    flag2 = _PAIRED | _READ2 | _REVERSE | (_PROPER if proper else 0)
    return [
        _make_read(
            name, ref_id, pos1, "50M", 50, flag=flag1,
            mate_ref_id=mate_ref, mate_start=pos2, nm=0,
        ),
        _make_read(
            name, mate_ref, pos2, "50M", 50, flag=flag2,
            mate_ref_id=ref_id, mate_start=pos1, nm=0,
        ),
    ]


def test_good_pairs_counted_once_from_read_one(tmp_path):
    reads = []
    for i in range(5):
        reads += _fr_pair(f"p{i}", 0, 10 + i * 20, 60 + i * 20)
    bam = _write_bam(tmp_path / "good.bam", reads)
    by_name = {c.name: c for c in compute_bam_metrics(bam, realistic_distance=500)}
    assert by_name["contigA"].good == 5
    assert by_name["contigA"].both_mapped == 5
    assert by_name["contigA"].fragments_mapped == 5
    assert by_name["contigA"].properpair == 5
    assert by_name["contigA"].bridges == 0


def test_mates_on_different_contigs_are_bridges(tmp_path):
    reads = _fr_pair("p0", 0, 10, 60, mate_ref_id=1)
    bam = _write_bam(tmp_path / "bridge.bam", reads)
    by_name = {c.name: c for c in compute_bam_metrics(bam, realistic_distance=500)}
    assert by_name["contigA"].bridges == 1
    assert by_name["contigA"].good == 0


def test_pairs_beyond_realistic_distance_are_not_good(tmp_path):
    reads = _fr_pair("p0", 0, 10, 350)
    bam = _write_bam(tmp_path / "far.bam", reads)
    by_name = {c.name: c for c in compute_bam_metrics(bam, realistic_distance=100)}
    assert by_name["contigA"].good == 0
    assert by_name["contigA"].both_mapped == 1


def test_wrong_orientation_is_not_good(tmp_path):
    # Both mates forward: neither FR nor RF.
    flag1 = _PAIRED | _READ1
    flag2 = _PAIRED | _READ2
    reads = [
        _make_read("p0", 0, 10, "50M", 50, flag=flag1, mate_ref_id=0,
                   mate_start=60, nm=0),
        _make_read("p0", 0, 60, "50M", 50, flag=flag2, mate_ref_id=0,
                   mate_start=10, nm=0),
    ]
    bam = _write_bam(tmp_path / "orient.bam", reads)
    by_name = {c.name: c for c in compute_bam_metrics(bam, realistic_distance=500)}
    assert by_name["contigA"].good == 0


def test_orphan_read_two_counts_as_a_mapped_fragment(tmp_path):
    read = _make_read(
        "orphan", 0, 10, "50M", 50,
        flag=_PAIRED | _READ2 | _MATE_UNMAPPED, nm=0,
    )
    bam = _write_bam(tmp_path / "orphan.bam", [read])
    by_name = {c.name: c for c in compute_bam_metrics(bam)}
    assert by_name["contigA"].fragments_mapped == 1
    assert by_name["contigA"].both_mapped == 0


# ---------------------------------------------------------------------------
# p_seq_true
# ---------------------------------------------------------------------------


def test_p_seq_true_rescaling(tmp_path):
    # scale = (100-35)/100 = 0.65; NM=0 -> (1 - .65)/.35 == 1.
    bam = _write_bam(
        tmp_path / "perfect.bam", [_make_read("r", 0, 10, "100M", 100, nm=0)]
    )
    by_name = {c.name: c for c in compute_bam_metrics(bam)}
    assert by_name["contigA"].p_seq_true() == pytest.approx(1.0)


def test_p_seq_true_is_zero_at_edit_distance_35(tmp_path):
    bam = _write_bam(
        tmp_path / "bad.bam", [_make_read("r", 0, 10, "100M", 100, nm=35)]
    )
    by_name = {c.name: c for c in compute_bam_metrics(bam)}
    assert by_name["contigA"].p_seq_true() == pytest.approx(0.0, abs=1e-12)


def test_p_seq_true_goes_negative_past_the_floor(tmp_path):
    """The C++ does not clamp; preserved so scores stay comparable."""
    bam = _write_bam(
        tmp_path / "worse.bam", [_make_read("r", 0, 10, "100M", 100, nm=70)]
    )
    by_name = {c.name: c for c in compute_bam_metrics(bam)}
    assert by_name["contigA"].p_seq_true() < 0


def test_contig_with_no_reads_reports_p_seq_true_of_one(tmp_path):
    bam = _write_bam(
        tmp_path / "empty.bam", [_make_read("r", 0, 10, "50M", 50, nm=0)]
    )
    by_name = {c.name: c for c in compute_bam_metrics(bam)}
    assert by_name["contigB"].reads_mapped == 0
    assert by_name["contigB"].p_seq_true() == 1.0


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------


def test_one_row_per_reference_in_header_order(tmp_path):
    bam = _write_bam(
        tmp_path / "order.bam", [_make_read("r", 1, 10, "50M", 50, nm=0)]
    )
    contigs = compute_bam_metrics(bam)
    assert [c.name for c in contigs] == ["contigA", "contigB"]
    assert [c.length for c in contigs] == [400, 300]


def test_row_keys_match_the_bam_read_columns(tmp_path):
    bam = _write_bam(
        tmp_path / "cols.bam", [_make_read("r", 0, 10, "50M", 50, nm=0)]
    )
    row = compute_bam_metrics(bam)[0].as_row()
    assert tuple(row.keys()) == CSV_COLUMNS


def test_fragment_estimate_is_positive_for_a_real_library(tmp_path):
    reads = []
    for i in range(50):
        reads += _fr_pair(f"p{i}", 0, 10 + i * 4, 160 + i * 4)
    # Read order, not coordinate order: the estimator pairs adjacent records.
    bam = _write_bam(tmp_path / "frag.bam", reads, sort=False)
    assert estimate_realistic_distance(bam) > 0


def test_fragment_estimate_needs_read_ordered_input(tmp_path):
    """Guards the ordering contract: a coordinate sort separates the mates."""
    reads = []
    for i in range(50):
        reads += _fr_pair(f"p{i}", 0, 10 + i * 4, 160 + i * 4)
    read_ordered = _write_bam(tmp_path / "ro.bam", reads, sort=False)
    coord_sorted = _write_bam(tmp_path / "cs.bam", reads, sort=True)
    assert estimate_realistic_distance(read_ordered) > 0
    assert estimate_realistic_distance(coord_sorted) == 0


def test_empty_contig_metrics_are_wellformed():
    contig = ContigMetrics(name="x", length=0)
    contig.calculate_uncovered_bases()
    contig.set_p_not_segmented()
    assert contig.bases_uncovered == 0
    assert 0.0 <= contig.p_not_segmented <= 1.0


def test_contig_metrics_compare_on_their_counters():
    """The private coverage fields are declared but stay out of __eq__.

    See SLOTTED_ACCUMULATOR: they became fields so the class could take
    __slots__, and a numpy array in __eq__ would raise rather than compare.
    """
    left = ContigMetrics(name="x", length=10)
    right = ContigMetrics(name="x", length=10)
    assert left == right
    left._counts[3] += 1
    assert left == right
    left.good += 1
    assert left != right


def test_contig_metrics_reject_an_unknown_attribute():
    contig = ContigMetrics(name="x", length=10)
    with pytest.raises(AttributeError):
        contig.reads_maped = 1  # noqa: B010 - the typo is the point


# ---------------------------------------------------------------------------
# Soft-clip accounting
#
# Reported during the run because the BAM it comes from is deleted by
# default, and it is what explains a drop in sCcov under snap 2.x.
# ---------------------------------------------------------------------------


def test_soft_clips_are_counted(tmp_path):
    reads = [
        _make_read("r1", 0, 50, "20S80M", 100, nm=0),
        _make_read("r2", 0, 60, "80M20S", 100, nm=0),
        _make_read("r3", 0, 70, "100M", 100, nm=0),
    ]
    bam = _write_bam(tmp_path / "clipstats.bam", reads)
    by_name = {c.name: c for c in compute_bam_metrics(bam)}
    contig = by_name["contigA"]
    assert contig.clipped_alignments == 2
    assert contig.clipped_bases == 40
    assert contig.leading_clipped_bases == 20   # only r1 clips at the start


def test_unclipped_alignments_count_zero(tmp_path):
    bam = _write_bam(
        tmp_path / "noclip2.bam", [_make_read("r", 0, 10, "50M10D40M", 90, nm=0)]
    )
    contig = {c.name: c for c in compute_bam_metrics(bam)}["contigA"]
    assert contig.clipped_alignments == 0
    assert contig.clipped_bases == 0


def test_clip_counting_does_not_disturb_coverage(tmp_path):
    """The counters must not move the reference cursor."""
    bam = _write_bam(
        tmp_path / "clipcov.bam", [_make_read("r", 0, 100, "25S75M", 100, nm=0)]
    )
    got = _coverage_by_name(compute_bam_metrics(bam))["contigA"]
    np.testing.assert_array_equal(got, _samtools_depth(bam)["contigA"])
    assert int(np.flatnonzero(got)[0]) == 100


# ---------------------------------------------------------------------------
# MALFORMED_RECORDS
# ---------------------------------------------------------------------------

_SEQ_NT16 = {base: code for code, base in enumerate("=ACMGRSVTWYHKDBN")}


def _raw_bam_record(name, ref_id, pos, cigar, seq, flag=0, mapq=60):
    """Encode one BAM record by hand.

    pysam validates on write, so a record with a CIGAR that disagrees with
    the read length -- the thing snap emits at a contig boundary -- cannot be
    produced through AlignedSegment.  See MALFORMED_RECORDS.
    """
    import struct

    qname = name.encode() + b"\0"
    packed_cigar = b"".join(struct.pack("<I", (n << 4) | op) for n, op in cigar)
    packed_seq = bytearray()
    for i in range(0, len(seq), 2):
        high = _SEQ_NT16[seq[i]]
        low = _SEQ_NT16[seq[i + 1]] if i + 1 < len(seq) else 0
        packed_seq.append(high << 4 | low)
    core = struct.pack(
        "<iiBBHHHiiii", ref_id, pos, len(qname), mapq, 4680,
        len(cigar), flag, len(seq), -1, -1, 0,
    )
    body = core + qname + packed_cigar + bytes(packed_seq) + b"\xff" * len(seq)
    return struct.pack("<i", len(body)) + body


def _write_raw_bam(path, records, refs=REFS):
    """Write a BGZF BAM from hand-encoded records."""
    import struct

    header = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in refs
    )
    out = bytearray(b"BAM\1")
    out += struct.pack("<i", len(header)) + header.encode()
    out += struct.pack("<i", len(refs))
    for name, length in refs:
        encoded = name.encode() + b"\0"
        out += struct.pack("<i", len(encoded)) + encoded + struct.pack("<i", length)
    for record in records:
        out += record

    plain = str(path) + ".plain"
    with open(plain, "wb") as handle:
        handle.write(bytes(out))
    pysam.tabix_compress(plain, str(path), force=True)
    return str(path)


def test_htslib_resumes_after_a_record_it_refuses(tmp_path):
    """The load-bearing claim behind skipping: htslib consumes the whole
    record before it checks the CIGAR, so the stream is already at the next
    record when it returns -4.  Asserted against real htslib, not a fake.
    """
    from pytransrate.bam_metrics import MalformedRecordStats, iter_alignments

    seq = "ACGT" * 25  # 100 bp
    bam_path = _write_raw_bam(
        tmp_path / "boundary.bam",
        [
            _raw_bam_record("good_1", 0, 10, [(100, 0)], seq),
            _raw_bam_record("broken", 0, 30, [(90, 0)], seq),      # 90M vs 100 bp
            _raw_bam_record("good_2", 0, 40, [(10, 4), (90, 0)], seq),
            _raw_bam_record("broken_2", 0, 60, [(47, 4), (48, 4)], seq),
            _raw_bam_record("good_3", 0, 70, [(100, 0)], seq),
        ],
    )

    stats = MalformedRecordStats()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        names = [read.query_name for read in iter_alignments(bam, bam_path, stats)]

    assert names == ["good_1", "good_2", "good_3"]  # nothing after is lost
    assert stats.skipped == 2


def test_skipped_records_are_reported_not_swallowed():
    """Skipped alignments are lost coverage, so the count must surface."""
    from pytransrate.bam_metrics import MalformedRecordStats

    stats = MalformedRecordStats(skipped=7)
    assert "7 alignment" in stats.describe()
    assert "unmapped" in stats.describe()


def test_a_run_of_unreadable_records_is_corruption_and_stops(tmp_path):
    """-4 also means truncation, where the stream is desynchronised and
    continuing would be nonsense.  Told apart by shape: snap's are isolated.
    """
    from pytransrate.bam_metrics import (
        MALFORMED_RUN_LIMIT,
        MalformedBamError,
        iter_alignments,
    )

    # A class, not a generator: pysam's iterator survives raising, which is
    # the whole reason skipping works.  A generator fake would not.
    class _Bam:
        calls = 0

        def fetch(self, **kwargs):
            return self

        def __iter__(self):
            return self

        def __next__(self):
            self.calls += 1
            if self.calls == 1:
                return "first record"
            raise OSError("error -4 while reading file")

    seen = []
    with pytest.raises(MalformedBamError) as caught:
        for read in iter_alignments(_Bam(), "aln.bam"):
            seen.append(read)

    assert seen == ["first record"]
    message = str(caught.value)
    assert "aln.bam" in message
    assert str(MALFORMED_RUN_LIMIT) in message
    assert "truncated or corrupt" in message


def test_the_run_counter_resets_on_a_good_record():
    """Scattered bad records must never add up to a corruption verdict."""
    from pytransrate.bam_metrics import (
        MALFORMED_RUN_LIMIT,
        MalformedRecordStats,
        iter_alignments,
    )

    class _Bam:
        index = -1

        def fetch(self, **kwargs):
            return self

        def __iter__(self):
            return self

        def __next__(self):
            self.index += 1
            if self.index >= MALFORMED_RUN_LIMIT * 3:
                raise StopIteration
            if self.index % 2:
                raise OSError("error -4 while reading file")
            return f"read_{self.index}"

    stats = MalformedRecordStats()
    reads = list(iter_alignments(_Bam(), "aln.bam", stats))

    assert len(reads) == MALFORMED_RUN_LIMIT * 3 // 2
    assert stats.skipped == MALFORMED_RUN_LIMIT * 3 // 2


def test_clean_iteration_is_unaffected(tmp_path):
    from pytransrate.bam_metrics import iter_alignments

    bam_path = _write_bam(
        tmp_path / "fine.bam", [_make_read("r", 0, 10, "50M", 50, nm=0)]
    )
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        assert len(list(iter_alignments(bam, bam_path))) == 1
