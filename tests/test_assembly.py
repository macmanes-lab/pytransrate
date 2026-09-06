"""Tests for FASTA loading, contig metrics, and assembly statistics."""

from __future__ import annotations

import pytest

from transrate.assembly import (
    BASIC_STATS_KEYS,
    Assembly,
    AssemblyError,
    parse_fasta,
)
from transrate.contig import SCORE_FLOOR, Contig


def _write(tmp_path, text, name="a.fa"):
    path = tmp_path / name
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def test_parses_multiline_sequences(tmp_path):
    path = _write(tmp_path, ">c1\nACGT\nACGT\n>c2\nTTTT\n")
    assert list(parse_fasta(path)) == [("c1", "ACGTACGT"), ("c2", "TTTT")]


def test_identifier_stops_at_whitespace_or_pipe(tmp_path):
    path = _write(tmp_path, ">c1 a description here\nACGT\n>c2|x|y\nACGT\n")
    assert [n for n, _ in parse_fasta(path)] == ["c1", "c2"]


def test_empty_sequence_is_rejected(tmp_path):
    path = _write(tmp_path, ">c1\n\n>c2\nACGT\n")
    with pytest.raises(AssemblyError, match="no sequence"):
        Assembly(path)


def test_duplicate_identifiers_are_rejected(tmp_path):
    # The Trinity '|' bug: two deflines collapsing to one identifier.
    path = _write(tmp_path, ">c1|a\nACGT\n>c1|b\nTTTT\n")
    with pytest.raises(AssemblyError, match="Non-unique"):
        Assembly(path)


def test_commas_in_names_are_rejected(tmp_path):
    path = _write(tmp_path, ">c1,x\nACGT\n")
    with pytest.raises(AssemblyError, match="commas"):
        Assembly(path)


def test_n_bases_counts_every_contig(tmp_path):
    path = _write(tmp_path, ">c1\n" + "A" * 100 + "\n>c2\n" + "C" * 50 + "\n")
    assert Assembly(path).n_bases == 150


# ---------------------------------------------------------------------------
# Contig scoring
# ---------------------------------------------------------------------------


def test_score_is_the_product_of_four_components():
    contig = Contig(name="c", seq="A" * 100)
    contig.set_uncovered_bases(0)
    contig.p_not_segmented = 0.5
    contig.p_good = 0.5
    contig.p_seq_true = 0.5
    assert contig.score == pytest.approx(1.0 * 0.5 * 0.5 * 0.5)


def test_score_components_are_floored():
    contig = Contig(name="c", seq="A" * 100)
    contig.set_uncovered_bases(100)  # p_bases_covered == 0
    contig.p_not_segmented = 0.0
    contig.p_good = 0.0
    contig.p_seq_true = 0.0
    assert contig.score == SCORE_FLOOR


def test_score_never_falls_below_the_floor():
    contig = Contig(name="c", seq="A" * 100)
    contig.set_uncovered_bases(0)
    contig.p_not_segmented = 0.02
    contig.p_good = 0.02
    contig.p_seq_true = 0.02
    assert contig.score == SCORE_FLOOR


def test_alt_scores_leave_one_component_out():
    contig = Contig(name="c", seq="A" * 100)
    contig.set_uncovered_bases(0)
    contig.p_not_segmented = 0.5
    contig.p_good = 0.4
    contig.p_seq_true = 0.2
    alt = contig.alt_scores()
    assert alt["cov"] == pytest.approx(0.5 * 0.4 * 0.2)
    assert alt["seg"] == pytest.approx(1.0 * 0.4 * 0.2)


def test_classification():
    contig = Contig(name="c", seq="A" * 100)
    contig.set_uncovered_bases(0)
    contig.p_not_segmented = contig.p_good = contig.p_seq_true = 1.0
    assert contig.classify(0.5) == "good"
    assert contig.classify(1.5) == "bad"


def test_trailing_semicolons_are_stripped_from_names():
    assert Contig(name="c1;", seq="ACGT").name == "c1"


# ---------------------------------------------------------------------------
# Assembly statistics -- see STATS_QUIRK
# ---------------------------------------------------------------------------


def test_basic_stats_key_order_is_stable(tmp_path):
    path = _write(tmp_path, ">c1\n" + "ACGT" * 100 + "\n")
    assert tuple(Assembly(path).basic_stats().keys()) == BASIC_STATS_KEYS


def test_nx_labels_track_ascending_walk(tmp_path):
    # Ten 1000bp contigs: every Nx is 1000.
    text = "".join(f">c{i}\n{'ACGT' * 250}\n" for i in range(10))
    stats = Assembly(_write(tmp_path, text)).basic_stats()
    for label in ("n90", "n70", "n50", "n30", "n10"):
        assert stats[label] == 1000, label


def test_nx_on_mixed_lengths(tmp_path):
    text = "".join(f">s{i}\n{'A' * 300}\n" for i in range(9))
    text += ">big\n" + "A" * 10000 + "\n"
    stats = Assembly(_write(tmp_path, text)).basic_stats()
    # Walking shortest-first, 10% of bases is reached among the short
    # contigs, so N90 is small and N10 is the long one.
    assert stats["n90"] == 300
    assert stats["n10"] == 10000
    assert stats["largest"] == 10000
    assert stats["smallest"] == 300


def test_short_contigs_counted_but_excluded_from_lengths(tmp_path):
    """STATS_QUIRK 1: <200bp contigs skew mean_len downward."""
    text = "".join(f">short{i}\n{'A' * 100}\n" for i in range(10))
    text += ">long\n" + "A" * 1000 + "\n"
    stats = Assembly(_write(tmp_path, text)).basic_stats()

    assert stats["n_seqs"] == 11
    assert stats["n_under_200"] == 10
    assert stats["n_bases"] == 2000  # all contigs
    # cumulative (1000, the long contig only) / 11 contigs
    assert stats["mean_len"] == pytest.approx(1000 / 11)


def test_nx_padding_uses_the_longest_contig(tmp_path):
    """STATS_QUIRK 2 and 3: unreached thresholds pad with the largest."""
    text = "".join(f">short{i}\n{'A' * 100}\n" for i in range(50))
    text += ">long\n" + "A" * 500 + "\n"
    stats = Assembly(_write(tmp_path, text)).basic_stats()
    # Only one contig clears 200bp, so at most one threshold is satisfied
    # per contig and the rest pad with the longest length.
    assert stats["n10"] == 500


def test_n_with_orf_threshold(tmp_path):
    # 150 codons of ATG with no stop -> orf_length 150 > 149.
    path = _write(tmp_path, ">c1\n" + "ATG" * 150 + "\n")
    assert Assembly(path).basic_stats()["n_with_orf"] == 1


def test_empty_assembly_stats_are_zeroed(tmp_path):
    path = _write(tmp_path, "")
    stats = Assembly(path).basic_stats()
    assert set(stats.keys()) == set(BASIC_STATS_KEYS)
    assert all(v == 0 for v in stats.values())


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def test_contig_metrics_gc_excludes_n_but_proportion_n_does_not(tmp_path):
    # 2 G, 2 C, 2 A, 2 T, 2 N -> gc = 4/8, proportion_n = 2/10.
    path = _write(tmp_path, ">c1\nGGCCAATTNN\n")
    metrics = Assembly(path).contig_metrics()
    assert metrics["gc"] == pytest.approx(0.5)
    assert metrics["bases_n"] == 2
    assert metrics["proportion_n"] == pytest.approx(0.2)


def test_contig_metrics_on_all_n(tmp_path):
    path = _write(tmp_path, ">c1\nNNNN\n")
    metrics = Assembly(path).contig_metrics()
    assert metrics["gc"] == 0.0
    assert metrics["proportion_n"] == 1.0


def test_good_contigs_counts_classified(tmp_path):
    text = "".join(f">c{i}\n{'ACGT' * 100}\n" for i in range(4))
    assembly = Assembly(_write(tmp_path, text))
    for i, (_, contig) in enumerate(assembly):
        contig.set_uncovered_bases(0)
        contig.p_not_segmented = contig.p_good = contig.p_seq_true = (
            1.0 if i < 3 else 0.0
        )
    assembly.classify_contigs(0.5)
    assert assembly.good_contigs == 3
