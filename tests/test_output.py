"""Pins the CSV column contract that ORP depends on positionally.

These are the highest-consequence assertions in the suite.  If ``score``
moves off index 8 of contigs.csv, ORP's orthofuser silently picks contigs by
the wrong column; if it moves off index 36 of assemblies.csv, the reported
TRANSRATE SCORE is silently some other statistic.  Neither failure raises.
"""

from __future__ import annotations

import csv

import pytest

from transrate.assembly import Assembly
from transrate.output import (
    assemblies_csv_columns,
    contigs_csv_columns,
    write_assemblies_csv,
    write_contigs_csv,
)


# ---------------------------------------------------------------------------
# The two indices ORP hard-codes
# ---------------------------------------------------------------------------


def test_contigs_csv_score_is_column_nine():
    """scripts/pick_best_contigs.py reads float(row[8])."""
    columns = contigs_csv_columns(with_reference=False, with_reads=True)
    assert columns[0] == "contig_name"
    assert columns[8] == "score"


def test_assemblies_csv_score_indices():
    """oyster.py reads rows[1][36] and rows[1][37]."""
    columns = assemblies_csv_columns(with_reference=False, with_reads=True)
    assert columns[0] == "assembly"
    assert columns[36] == "score"
    assert columns[37] == "optimal_score"


def test_assemblies_csv_block_layout():
    columns = assemblies_csv_columns(with_reference=False, with_reads=True)
    # 1 assembly + 15 basic + 3 composition + 17 read + 4 trailing
    assert len(columns) == 40
    assert columns[1] == "n_seqs"
    assert columns[15] == "n10"
    assert columns[16] == "gc"
    assert columns[19] == "fragments"
    assert columns[35] == "p_contigs_segmented"
    assert columns[-1] == "weighted"


def test_contigs_csv_full_layout():
    assert contigs_csv_columns(with_reference=False, with_reads=True) == (
        "contig_name",
        "length",
        "prop_gc",
        "orf_length",
        "in_bridges",
        "p_good",
        "p_bases_covered",
        "p_seq_true",
        "score",
        "p_not_segmented",
        "eff_length",
        "eff_count",
        "tpm",
        "coverage",
        "sCnuc",
        "sCcov",
        "sCord",
        "sCseg",
    )


def test_reference_columns_precede_read_columns():
    columns = contigs_csv_columns(with_reference=True, with_reads=True)
    assert columns.index("has_crb") < columns.index("in_bridges")
    # Adding reference columns must not disturb the leading block.
    assert columns[:4] == ("contig_name", "length", "prop_gc", "orf_length")


def test_contig_only_run_omits_read_columns():
    columns = contigs_csv_columns(with_reference=False, with_reads=False)
    assert columns == ("contig_name", "length", "prop_gc", "orf_length")


# ---------------------------------------------------------------------------
# Round-tripping actual files
# ---------------------------------------------------------------------------


@pytest.fixture
def assembly(tmp_path):
    fasta = tmp_path / "a.fa"
    fasta.write_text(
        ">c1 some description\n" + "ATGAAACCCGGGTTT" * 20 + "\n"
        ">c2|withpipe\n" + "ACGT" * 100 + "\n"
    )
    return Assembly(fasta)


def test_contigs_csv_round_trip(tmp_path, assembly):
    for _, contig in assembly:
        contig.p_good = 0.9
        contig.p_seq_true = 0.95
        contig.p_not_segmented = 0.99
        contig.set_uncovered_bases(0)
        contig.tpm = 12.5

    path = tmp_path / "contigs.csv"
    write_contigs_csv(assembly, str(path))

    rows = list(csv.reader(open(path)))
    assert rows[0] == list(contigs_csv_columns())
    assert len(rows) == 3
    # The score column must parse as a float, the way ORP reads it.
    for row in rows[1:]:
        assert 0.0 < float(row[8]) <= 1.0


def test_contig_names_are_split_on_space_and_pipe(assembly):
    assert [name for name, _ in assembly] == ["c1", "c2"]


def test_assemblies_csv_round_trip(tmp_path):
    result = {c: 0 for c in assemblies_csv_columns()}
    result["assembly"] = "/path/to/a.fa"
    result["score"] = 0.515583333
    result["optimal_score"] = 0.526854321

    path = tmp_path / "assemblies.csv"
    write_assemblies_csv([result], str(path))

    rows = list(csv.reader(open(path)))
    assert rows[0][36] == "score"
    assert float(rows[1][36]) == pytest.approx(0.51558, abs=1e-5)
    assert float(rows[1][37]) == pytest.approx(0.52685, abs=1e-5)


def test_missing_metrics_render_as_empty_not_crash(tmp_path):
    path = tmp_path / "sparse.csv"
    write_assemblies_csv([{"assembly": "x.fa"}], str(path))
    rows = list(csv.reader(open(path)))
    assert rows[1][0] == "x.fa"
    assert rows[1][36] == ""
