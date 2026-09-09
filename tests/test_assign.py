"""Tests for in-memory fragment assignment.

This replaces salmon 0.8.2's --sampleOut. It is a redesign rather than a
port, so these tests pin the properties we actually want -- determinism,
preference for proper pairs, and sensitivity to both alignment quality and
abundance -- rather than agreement with a reference implementation that no
longer exists.
"""

from __future__ import annotations

import pysam
import pytest

from pytransrate.assign import (
    ORPHAN_EDIT_FRACTION,
    assign_fragments,
    group_by_fragment,
    score_candidates,
)
from pytransrate.bam_metrics import decode

REFS = ["txA", "txB"]

_PAIRED = 1
_REVERSE = 16
_MATE_REVERSE = 32
_READ1 = 64
_READ2 = 128
_SECONDARY = 256


def _aln(name, ref_id, start, *, nm=0, read1=True, secondary=False, length=100):
    read = pysam.AlignedSegment()
    read.query_name = name
    flag = _PAIRED | (_READ1 | _MATE_REVERSE if read1 else _READ2 | _REVERSE)
    if secondary:
        flag |= _SECONDARY
    read.flag = flag
    read.reference_id = ref_id
    read.reference_start = start
    read.mapping_quality = 60
    read.cigarstring = f"{length}M"
    read.query_sequence = "A" * length
    read.query_qualities = pysam.qualitystring_to_array("I" * length)
    read.set_tag("NM", nm, value_type="i")
    return read


def _decoded(reads):
    """A batch as score_candidates now takes it. See DECODE_ONCE."""
    return [decode(read) for read in reads]


def _expr(**counts):
    return {name: {"eff_count": c, "eff_len": 1000, "tpm": c}
            for name, c in counts.items()}


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------


def test_group_by_fragment_batches_contiguous_records():
    reads = [
        _aln("f1", 0, 10), _aln("f1", 0, 200, read1=False),
        _aln("f2", 1, 10), _aln("f2", 1, 200, read1=False),
    ]
    groups = list(group_by_fragment(reads))
    assert [name for name, _ in groups] == ["f1", "f2"]
    assert [len(batch) for _, batch in groups] == [2, 2]


def test_group_by_fragment_yields_the_final_batch():
    groups = list(group_by_fragment([_aln("only", 0, 10)]))
    assert len(groups) == 1


def test_group_by_fragment_on_empty_input():
    assert list(group_by_fragment([])) == []


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def test_better_alignment_wins_at_equal_abundance():
    batch = [
        _aln("f", 0, 10, nm=0), _aln("f", 0, 200, nm=0, read1=False),
        _aln("f", 1, 10, nm=5, secondary=True),
        _aln("f", 1, 200, nm=5, read1=False, secondary=True),
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 10.0, "txB": 10.0})
    assert scored[0][0] > scored[1][0]


def test_higher_abundance_wins_at_equal_alignment_quality():
    batch = [
        _aln("f", 0, 10, nm=2), _aln("f", 0, 200, nm=2, read1=False),
        _aln("f", 1, 10, nm=2, secondary=True),
        _aln("f", 1, 200, nm=2, read1=False, secondary=True),
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 1.0, "txB": 500.0})
    assert scored[1][0] > scored[0][0]


def test_abundance_can_be_outweighed_by_a_much_better_alignment():
    batch = [
        _aln("f", 0, 10, nm=0), _aln("f", 0, 200, nm=0, read1=False),
        _aln("f", 1, 10, nm=30, secondary=True),
        _aln("f", 1, 200, nm=30, read1=False, secondary=True),
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 1.0, "txB": 100.0})
    assert scored[0][0] > scored[1][0]


def test_proper_pair_beats_an_orphan_placement():
    batch = [
        _aln("f", 0, 10, nm=1), _aln("f", 0, 200, nm=1, read1=False),
        _aln("f", 1, 10, nm=0, secondary=True),  # better, but one mate only
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 10.0, "txB": 10.0})
    assert scored[0][0] > scored[1][0]


@pytest.mark.parametrize("nm", [0, 1, 2, 3, 4, 5, 6, 7])
def test_pair_beats_a_perfect_orphan_up_to_the_crossover(nm):
    """A reasonably-aligned pair outranks a single perfect mate."""
    batch = [
        _aln("f", 0, 10, nm=nm), _aln("f", 0, 200, nm=nm, read1=False),
        _aln("f", 1, 10, nm=0, secondary=True),  # perfect, but one mate only
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 10.0, "txB": 10.0})
    assert scored[0][0] > scored[1][0]


@pytest.mark.parametrize("nm", [8, 12, 20])
def test_badly_diverged_pair_loses_to_a_perfect_orphan(nm):
    """Deliberate: past ~7.5% per-mate divergence the pair is not evidence
    the fragment came from that transcript. See ORPHAN_EDIT_FRACTION."""
    batch = [
        _aln("f", 0, 10, nm=nm), _aln("f", 0, 200, nm=nm, read1=False),
        _aln("f", 1, 10, nm=0, secondary=True),
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 10.0, "txB": 10.0})
    assert scored[1][0] > scored[0][0]


def test_orphan_charge_cancels_when_no_candidate_has_the_mate():
    """A genuinely unpaired fragment ranks on the mate it does have."""
    batch = [
        _aln("f", 0, 10, nm=4),
        _aln("f", 1, 10, nm=0, secondary=True),
    ]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 10.0, "txB": 10.0})
    assert scored[1][0] > scored[0][0]


def test_orphan_edit_fraction_is_a_fraction():
    assert 0.0 < ORPHAN_EDIT_FRACTION < 1.0


@pytest.mark.parametrize(
    "settings",
    [
        {"orphan_edit_fraction": 0.05},
        {"orphan_edit_fraction": 0.40},
        {"error_rate": 0.001},
        {"error_rate": 0.10},
    ],
)
def test_the_orphan_charge_follows_its_settings(settings):
    """The charge is memoised, so it has to be keyed on what it depends on."""
    batch = [
        _aln("f", 0, 10, nm=0), _aln("f", 0, 200, nm=0, read1=False),
        _aln("f", 1, 10, nm=0, secondary=True),  # one mate only: charged
    ]
    decoded = _decoded(batch)
    default = score_candidates(decoded, REFS, {"txA": 10.0, "txB": 10.0})
    changed = score_candidates(decoded, REFS, {"txA": 10.0, "txB": 10.0}, **settings)

    # txA explains both mates, so only txB's score moves with the charge.
    assert changed[1][0] != default[1][0]


def test_the_orphan_charge_is_not_levied_on_a_complete_candidate():
    """A candidate explaining every mate is scored on its alignments alone."""
    batch = [_aln("f", 0, 10, nm=3), _aln("f", 0, 200, nm=3, read1=False)]
    decoded = _decoded(batch)
    scored = score_candidates(decoded, REFS, {"txA": 10.0})
    loosened = score_candidates(
        decoded, REFS, {"txA": 10.0}, orphan_edit_fraction=0.9
    )
    assert scored[0][0] == loosened[0][0]


def test_unmapped_records_are_ignored():
    read = _aln("f", 0, 10)
    read.flag = read.flag | 4
    read.reference_id = -1
    assert score_candidates(_decoded([read]), REFS, {}) == {}


def test_transcripts_absent_from_salmon_stay_reachable():
    """The pseudocount keeps zero-abundance transcripts scoreable."""
    batch = [_aln("f", 1, 10, nm=0), _aln("f", 1, 200, nm=0, read1=False)]
    scored = score_candidates(_decoded(batch), REFS, {"txA": 99.0})
    assert 1 in scored


# ---------------------------------------------------------------------------
# Assignment
# ---------------------------------------------------------------------------


def test_one_transcript_survives_per_fragment():
    reads = [
        _aln("f", 0, 10, nm=0), _aln("f", 0, 200, nm=0, read1=False),
        _aln("f", 1, 10, nm=5, secondary=True),
        _aln("f", 1, 200, nm=5, read1=False, secondary=True),
    ]
    out = list(assign_fragments(reads, REFS, _expr(txA=10.0, txB=10.0)))
    assert len(out) == 2
    assert {r.reference_id for r in out} == {0}


def test_secondary_flag_is_cleared_on_survivors():
    reads = [
        _aln("f", 0, 10, nm=9), _aln("f", 0, 200, nm=9, read1=False),
        _aln("f", 1, 10, nm=0, secondary=True),
        _aln("f", 1, 200, nm=0, read1=False, secondary=True),
    ]
    out = list(assign_fragments(reads, REFS, _expr(txA=1.0, txB=1.0)))
    assert {r.reference_id for r in out} == {1}
    assert not any(r.is_secondary for r in out)


def test_uniquely_mapping_fragments_pass_through():
    reads = [_aln("f", 0, 10), _aln("f", 0, 200, read1=False)]
    out = list(assign_fragments(reads, REFS, _expr(txA=10.0)))
    assert len(out) == 2


def test_assignment_is_deterministic_regardless_of_record_order():
    a = _aln("f", 0, 10, nm=2)
    b = _aln("f", 0, 200, nm=2, read1=False)
    c = _aln("f", 1, 10, nm=2, secondary=True)
    d = _aln("f", 1, 200, nm=2, read1=False, secondary=True)
    expr = _expr(txA=50.0, txB=50.0)

    first = [r.reference_id for r in assign_fragments([a, b, c, d], REFS, expr)]
    # Rebuild in the opposite candidate order; same fragment either way.
    a2 = _aln("f", 1, 10, nm=2)
    b2 = _aln("f", 1, 200, nm=2, read1=False)
    c2 = _aln("f", 0, 10, nm=2, secondary=True)
    d2 = _aln("f", 0, 200, nm=2, read1=False, secondary=True)
    second = [r.reference_id for r in assign_fragments([a2, b2, c2, d2], REFS, expr)]

    # Exact tie on score and prior -> the lexicographically smaller name wins,
    # not whichever the aligner happened to list first.
    assert set(first) == set(second) == {0}


def test_assignment_without_expression_uses_alignment_quality_only():
    reads = [
        _aln("f", 0, 10, nm=7), _aln("f", 0, 200, nm=7, read1=False),
        _aln("f", 1, 10, nm=0, secondary=True),
        _aln("f", 1, 200, nm=0, read1=False, secondary=True),
    ]
    out = list(assign_fragments(reads, REFS, expression=None))
    assert {r.reference_id for r in out} == {1}


def test_multiple_fragments_are_assigned_independently():
    reads = [
        _aln("f1", 0, 10, nm=0), _aln("f1", 0, 200, nm=0, read1=False),
        _aln("f1", 1, 10, nm=5, secondary=True),
        _aln("f1", 1, 200, nm=5, read1=False, secondary=True),
        _aln("f2", 0, 10, nm=5), _aln("f2", 0, 200, nm=5, read1=False),
        _aln("f2", 1, 10, nm=0, secondary=True),
        _aln("f2", 1, 200, nm=0, read1=False, secondary=True),
    ]
    out = list(assign_fragments(reads, REFS, _expr(txA=10.0, txB=10.0)))
    by_frag = {}
    for read in out:
        by_frag.setdefault(read.query_name, set()).add(read.reference_id)
    assert by_frag == {"f1": {0}, "f2": {1}}
