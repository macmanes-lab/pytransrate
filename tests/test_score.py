"""Tests for assembly scoring and cutoff optimisation."""

from __future__ import annotations

import csv
import math

import pytest

from pytransrate.assembly import Assembly
from pytransrate.score import ScoreOptimiser, geomean


def _assembly(tmp_path, specs):
    """Build an assembly whose contigs have the given (p_good, good, tpm)."""
    text = "".join(f">c{i}\n{'ACGT' * 100}\n" for i in range(len(specs)))
    path = tmp_path / "a.fa"
    path.write_text(text)
    assembly = Assembly(path)
    for (name, contig), (p_good, good, tpm) in zip(assembly, specs):
        contig.set_uncovered_bases(0)
        contig.p_not_segmented = 1.0
        contig.p_seq_true = 1.0
        contig.p_good = p_good
        contig.good = good
        contig.tpm = tpm
    return assembly


# ---------------------------------------------------------------------------
# geomean
# ---------------------------------------------------------------------------


def test_geomean_of_identical_values():
    assert geomean([0.5] * 10) == pytest.approx(0.5)


def test_geomean_is_not_the_arithmetic_mean():
    values = [0.1, 1.0]
    assert geomean(values) == pytest.approx(math.sqrt(0.1))
    assert geomean(values) < sum(values) / len(values)


def test_geomean_of_empty_is_zero():
    assert geomean([]) == 0.0


# ---------------------------------------------------------------------------
# raw and weighted scores
# ---------------------------------------------------------------------------


def test_raw_score_scales_geomean_by_good_mapping_rate(tmp_path):
    assembly = _assembly(tmp_path, [(0.5, 10, 1.0)] * 4)
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=40)
    # every contig scores 0.5; good rate is 0.4
    assert optimiser.raw_score() == pytest.approx(0.5 * 0.4)


def test_weighted_score_uses_expression(tmp_path):
    assembly = _assembly(tmp_path, [(0.5, 10, 2.0), (0.5, 10, 0.0)])
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=50)
    # mean of (0.5*2.0, 0.5*0.0) == 0.5, times the 0.5 good rate
    assert optimiser.weighted_score() == pytest.approx(0.5 * 0.5)


def test_scores_are_zero_without_fragments(tmp_path):
    assembly = _assembly(tmp_path, [(0.5, 10, 1.0)])
    optimiser = ScoreOptimiser(assembly=assembly, fragments=0, good=0)
    assert optimiser.raw_score() == 0.0
    assert optimiser.weighted_score() == 0.0
    assert optimiser.optimal_score() == (0.0, 0.0)


def test_empty_assembly_scores_zero(tmp_path):
    path = tmp_path / "empty.fa"
    path.write_text("")
    optimiser = ScoreOptimiser(
        assembly=Assembly(path), fragments=100, good=50
    )
    assert optimiser.raw_score() == 0.0
    assert optimiser.optimal_score() == (0.0, 0.0)


# ---------------------------------------------------------------------------
# cutoff optimisation
# ---------------------------------------------------------------------------


def test_optimal_cutoff_discards_the_bad_contigs(tmp_path):
    # Three strong contigs plus three weak ones contributing no good reads.
    specs = [(1.0, 30, 1.0)] * 3 + [(0.01, 0, 1.0)] * 3
    assembly = _assembly(tmp_path, specs)
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=90)

    optimal, cutoff = optimiser.optimal_score()
    # Cutting the weak contigs raises the geometric mean without losing
    # good reads, so the optimum should beat the uncut score.
    assert optimal > optimiser.raw_score()
    assert 0.0 < cutoff <= 1.0


def test_optimal_score_is_cached(tmp_path):
    assembly = _assembly(tmp_path, [(1.0, 10, 1.0)] * 3)
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=30)
    first = optimiser.optimal_score()
    assert optimiser.optimal_score() is not None
    assert optimiser.optimal_score() == first


def test_optimisation_csv_is_written(tmp_path):
    assembly = _assembly(tmp_path, [(1.0, 10, 1.0), (0.5, 5, 1.0)])
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=15)
    path = tmp_path / "opt.csv"
    optimiser.optimal_score(csv_path=str(path))

    rows = list(csv.reader(open(path)))
    assert rows[0] == ["cutoff", "assembly_score"]
    assert len(rows) >= 2
    for row in rows[1:]:
        float(row[0])
        float(row[1])


def test_contigs_sharing_a_score_collapse_to_one_cutoff(tmp_path):
    """The Ruby keyed cutoffs by contig score, so ties collapse."""
    assembly = _assembly(tmp_path, [(0.5, 10, 1.0)] * 5)
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=50)
    path = tmp_path / "opt.csv"
    optimiser.optimal_score(csv_path=str(path))
    rows = list(csv.reader(open(path)))
    assert len(rows) == 2  # header plus one distinct cutoff


def test_final_contig_does_not_divide_by_zero(tmp_path):
    """The Ruby produced Inf/NaN on the last contig; we skip it instead."""
    assembly = _assembly(tmp_path, [(1.0, 10, 1.0)] * 3)
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=30)
    optimal, cutoff = optimiser.optimal_score()
    assert math.isfinite(optimal)
    assert math.isfinite(cutoff)


def test_single_contig_assembly(tmp_path):
    assembly = _assembly(tmp_path, [(1.0, 10, 1.0)])
    optimiser = ScoreOptimiser(assembly=assembly, fragments=100, good=10)
    optimal, cutoff = optimiser.optimal_score()
    # The only contig is the last one, so there is no valid cut.
    assert optimal == 0.0
    assert cutoff == 0.0
    assert optimiser.raw_score() == pytest.approx(0.1)
