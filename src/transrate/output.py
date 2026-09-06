"""CSV output, and the column contract that goes with it.

Port of ``Cmdline#write_contig_csv`` and ``Cmdline#write_assembly_csv``.

The column *order* here is an interface, not a presentation detail.  ORP
reads both files positionally:

* ``scripts/pick_best_contigs.py`` takes the contig id from column 1 and the
  score from column 9 (index 8) of ``contigs.csv``.
* ``oyster.py`` takes the transrate score from index 36 and the optimal score
  from index 37 of ``assemblies.csv``.

:func:`contigs_csv_columns` and :func:`assemblies_csv_columns` derive those
orders from the same key lists the metrics modules use, and
``tests/test_output.py`` pins the indices.  Change either at your peril.
"""

from __future__ import annotations

import csv
from collections import OrderedDict

from transrate.assembly import BASIC_STATS_KEYS, CONTIG_METRICS_KEYS

__all__ = [
    "READ_STATS_KEYS",
    "CONTIG_BASIC_KEYS",
    "CONTIG_READ_KEYS",
    "CONTIG_COMPARATIVE_KEYS",
    "ASSEMBLY_TRAILING_KEYS",
    "contigs_csv_columns",
    "assemblies_csv_columns",
    "write_contigs_csv",
    "write_assemblies_csv",
]

#: ``ReadMetrics#read_stats`` key order.
READ_STATS_KEYS = (
    "fragments",
    "fragments_mapped",
    "p_fragments_mapped",
    "good_mappings",
    "p_good_mapping",
    "bad_mappings",
    "potential_bridges",
    "bases_uncovered",
    "p_bases_uncovered",
    "contigs_uncovbase",
    "p_contigs_uncovbase",
    "contigs_uncovered",
    "p_contigs_uncovered",
    "contigs_lowcovered",
    "p_contigs_lowcovered",
    "contigs_segmented",
    "p_contigs_segmented",
)

#: ``Contig#basic_metrics``, prefixed with the name column.
CONTIG_BASIC_KEYS = ("contig_name", "length", "prop_gc", "orf_length")

#: ``Contig#read_metrics``.
CONTIG_READ_KEYS = (
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

#: ``Contig#comparative_metrics``, present only with ``--reference``.
CONTIG_COMPARATIVE_KEYS = ("has_crb", "reference_coverage", "hits")

#: Appended to every assemblies.csv row after the metric blocks.
ASSEMBLY_TRAILING_KEYS = ("score", "optimal_score", "cutoff", "weighted")


def contigs_csv_columns(
    with_reference: bool = False, with_reads: bool = True
) -> tuple[str, ...]:
    """Column order for ``contigs.csv``.

    Reference columns come before read columns, matching the order the Ruby
    merged them in ``write_contig_csv``.
    """
    columns = list(CONTIG_BASIC_KEYS)
    if with_reference:
        columns += list(CONTIG_COMPARATIVE_KEYS)
    if with_reads:
        columns += list(CONTIG_READ_KEYS)
    return tuple(columns)


def assemblies_csv_columns(
    with_reference: bool = False, with_reads: bool = True
) -> tuple[str, ...]:
    """Column order for ``assemblies.csv``.

    ``assembly`` is pulled to the front; everything else follows the order in
    which ``analyse_assembly`` merged the metric hashes.
    """
    columns = ["assembly"]
    columns += list(BASIC_STATS_KEYS)
    columns += list(CONTIG_METRICS_KEYS)
    if with_reads:
        columns += list(READ_STATS_KEYS)
    if with_reference:
        columns += list(_COMPARATIVE_STATS_KEYS)
    columns += list(ASSEMBLY_TRAILING_KEYS)
    return tuple(columns)


#: ``ComparativeMetrics#comp_stats`` key order. Unused without --reference.
_COMPARATIVE_STATS_KEYS = (
    "CRBB_hits",
    "n_contigs_with_CRBB",
    "p_contigs_with_CRBB",
    "rbh_per_reference",
    "n_refs_with_CRBB",
    "p_refs_with_CRBB",
    "cov25",
    "p_cov25",
    "cov50",
    "p_cov50",
    "cov75",
    "p_cov75",
    "cov85",
    "p_cov85",
    "cov95",
    "p_cov95",
    "reference_coverage",
)


def _round_floats(value, places: int):
    if isinstance(value, float):
        return round(value, places)
    return value


def write_contigs_csv(
    assembly,
    path: str,
    with_reference: bool = False,
    with_reads: bool = True,
) -> None:
    """Write per-contig metrics. Floats are rounded to 6 places, as in Ruby."""
    columns = contigs_csv_columns(with_reference, with_reads)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for name, contig in assembly:
            row = OrderedDict()
            row["contig_name"] = name
            row.update(contig.basic_metrics())
            if with_reference:
                row.update(contig.comparative_metrics())
            if with_reads:
                row.update(contig.read_metrics())
            writer.writerow([_round_floats(row[c], 6) for c in columns])


def write_assemblies_csv(
    results,
    path: str,
    with_reference: bool = False,
    with_reads: bool = True,
) -> None:
    """Write one row per assembly. Floats are rounded to 5 places, as in Ruby."""
    columns = assemblies_csv_columns(with_reference, with_reads)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for result in results:
            writer.writerow(
                [_round_floats(result.get(c, ""), 5) for c in columns]
            )
