#!/usr/bin/env python3
"""Compare transrate runs of one assembly under different mapping settings.

Answers two separate questions that are easy to conflate:

1. Does the setting change the score at all, and by how much?
2. Is that change bigger than the noise floor?

Question 2 is the one that decides the answer, and it needs a *replicate*:
two runs at the same settings on the same data.  snap and salmon are both
multithreaded and neither promises bit-identical output run to run
(amplab/snap#72 is literally "Inconsistent BAM output from identical
inputs"), so a delta smaller than the replicate delta says nothing.  Pass
``--replicate`` and every reported difference is scaled against it.

Usage::

    compare_transrate_runs.py --baseline default \\
        default=runs/orp_default \\
        replicate=runs/orp_default_rep \\
        mpc0=runs/orp_mpc0 \\
        mpc2=runs/orp_mpc2 \\
        --replicate replicate

Each directory is a transrate ``-o`` output directory holding
``assemblies.csv`` and ``contigs.csv``.  One invocation per assembly.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

# Columns computed from the FASTA alone.  They cannot move when only the
# mapping flags change, so a difference here means the runs were not given
# the same assembly and nothing below is worth reading.
try:  # keep in sync with the package when it is importable
    from pytransrate.assembly import BASIC_STATS_KEYS, CONTIG_METRICS_KEYS

    ASSEMBLY_SEQUENCE_KEYS = ("assembly",) + BASIC_STATS_KEYS + CONTIG_METRICS_KEYS
except ImportError:  # running somewhere pytransrate is not installed
    ASSEMBLY_SEQUENCE_KEYS = (
        "assembly", "n_seqs", "smallest", "largest", "n_bases", "mean_len",
        "n_under_200", "n_over_1k", "n_over_10k", "n_with_orf",
        "mean_orf_percent", "n90", "n70", "n50", "n30", "n10",
        "gc", "bases_n", "proportion_n",
    )

CONTIG_SEQUENCE_KEYS = ("length", "prop_gc", "orf_length")

#: Reported per run.  The first three are what ORP consumes; the rest are
#: there to explain any move in them.
ASSEMBLY_METRICS = (
    "score",
    "optimal_score",
    "cutoff",
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
    "contigs_uncovered",
    "contigs_lowcovered",
    "contigs_segmented",
)

#: score is the product of these four, each floored at 0.01 (contig.py).
#: Decomposing a score delta across them says *which* term moved.
SCORE_COMPONENTS = ("sCnuc", "sCcov", "sCord", "sCseg")

#: Below this a difference is float noise in a 6-decimal CSV, not a change.
TOL = 1e-6


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def load_run(label: str, directory: Path) -> dict:
    assemblies = directory / "assemblies.csv"
    contigs = directory / "contigs.csv"
    for path in (assemblies, contigs):
        if not path.exists():
            sys.exit(f"{label}: missing {path}")

    rows = read_csv(assemblies)
    if len(rows) != 1:
        sys.exit(
            f"{label}: {assemblies} has {len(rows)} rows; this script compares "
            "one assembly at a time, so run transrate once per assembly"
        )
    return {
        "label": label,
        "dir": directory,
        "assembly": rows[0],
        "contigs": {row["contig_name"]: row for row in read_csv(contigs)},
    }


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def column(run: dict, names: list[str], key: str) -> np.ndarray:
    """One contigs.csv column as floats, ordered by ``names``."""
    return np.array([as_float(run["contigs"][name][key]) for name in names])


# -- 1. validity ------------------------------------------------------------


def check_validity(runs: list[dict], baseline: dict) -> bool:
    """The runs must differ only in mapping. Anything else invalidates them."""
    print("=" * 74)
    print("1. VALIDITY  (sequence-only outputs must be identical)")
    print("=" * 74)

    ok = True
    base_names = set(baseline["contigs"])
    for run in runs:
        if run is baseline:
            continue
        problems = []

        for key in ASSEMBLY_SEQUENCE_KEYS:
            if key == "assembly":
                continue  # the path differs by design
            if key in baseline["assembly"] and (
                run["assembly"].get(key) != baseline["assembly"].get(key)
            ):
                problems.append(
                    f"assemblies.csv {key}: "
                    f"{baseline['assembly'].get(key)} -> {run['assembly'].get(key)}"
                )

        names = set(run["contigs"])
        if names != base_names:
            problems.append(
                f"contig set differs: {len(base_names - names)} missing, "
                f"{len(names - base_names)} extra"
            )
        else:
            for key in CONTIG_SEQUENCE_KEYS:
                differing = sum(
                    1
                    for name in base_names
                    if run["contigs"][name][key] != baseline["contigs"][name][key]
                )
                if differing:
                    problems.append(f"contigs.csv {key}: {differing} contigs differ")

        if problems:
            ok = False
            print(f"  FAIL  {run['label']}")
            for problem in problems:
                print(f"          {problem}")
        else:
            print(f"  ok    {run['label']}")

    if not ok:
        print(
            "\n  These runs did not analyse the same assembly. Fix that before\n"
            "  reading anything below -- the score differences are not the\n"
            "  mapping flags."
        )
    print()
    return ok


# -- 2. assembly level ------------------------------------------------------


def report_assembly_level(runs: list[dict], baseline: dict, replicate: dict | None):
    print("=" * 74)
    print("2. ASSEMBLY-LEVEL METRICS  (delta vs baseline)")
    print("=" * 74)

    others = [run for run in runs if run is not baseline]
    width = max([len(run["label"]) for run in runs] + [22]) + 2

    header = f"  {'metric':<22}{'baseline':>14}"
    for run in others:
        header += f"{run['label']:>{width}}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for key in ASSEMBLY_METRICS:
        if key not in baseline["assembly"]:
            continue
        base = as_float(baseline["assembly"][key])
        line = f"  {key:<22}{base:>14.6g}"
        for run in others:
            value = as_float(run["assembly"][key])
            delta = value - base
            line += f"{f'{value:.6g} ({delta:+.3g})':>{width}}"
        print(line)

    if replicate is not None:
        print()
        print("  Noise floor from the replicate, and each run's delta as a")
        print("  multiple of it. Under ~1x means indistinguishable from rerunning.")
        print()
        print(f"  {'metric':<22}{'|replicate delta|':>20}", end="")
        for run in others:
            if run is replicate:
                continue
            print(f"{run['label']:>14}", end="")
        print()
        print("  " + "-" * 72)
        for key in ASSEMBLY_METRICS:
            if key not in baseline["assembly"]:
                continue
            base = as_float(baseline["assembly"][key])
            noise = abs(as_float(replicate["assembly"][key]) - base)
            line = f"  {key:<22}{noise:>20.4g}"
            for run in others:
                if run is replicate:
                    continue
                delta = abs(as_float(run["assembly"][key]) - base)
                if noise > 0:
                    line += f"{f'{delta / noise:.1f}x':>14}"
                elif delta > 0:
                    line += f"{'inf':>14}"
                else:
                    line += f"{'0':>14}"
            print(line)
    else:
        print()
        print("  No --replicate given, so there is no noise floor and no way to")
        print("  tell a real effect from run-to-run variation. Add one.")
    print()


# -- 3. per contig ----------------------------------------------------------


def paired_summary(base: np.ndarray, other: np.ndarray) -> dict:
    delta = other - base
    changed = np.abs(delta) > TOL
    up = int(np.sum(delta > TOL))
    down = int(np.sum(delta < -TOL))
    moved = up + down
    return {
        "n": delta.size,
        "changed": int(changed.sum()),
        "p_changed": changed.sum() / delta.size if delta.size else 0.0,
        "mean": float(np.mean(delta)),
        "se": float(np.std(delta, ddof=1) / math.sqrt(delta.size))
        if delta.size > 1
        else 0.0,
        "median_changed": float(np.median(delta[changed])) if changed.any() else 0.0,
        "max_abs": float(np.max(np.abs(delta))) if delta.size else 0.0,
        "up": up,
        "down": down,
        # Sign test on the contigs that moved. With 100k contigs this goes
        # significant on a trivial shift, so read it as direction only --
        # the size of the effect is mean/median, and whether it matters is
        # the replicate comparison.
        "z": (up - down) / math.sqrt(moved) if moved else 0.0,
    }


def print_summary(label: str, summary: dict):
    lo = summary["mean"] - 1.96 * summary["se"]
    hi = summary["mean"] + 1.96 * summary["se"]
    print(
        f"    {label:<12}"
        f"changed {summary['changed']:>7} / {summary['n']:<7}"
        f"({summary['p_changed']:6.2%})  "
        f"mean {summary['mean']:+.3e} [{lo:+.2e},{hi:+.2e}]  "
        f"median|chg {summary['median_changed']:+.3e}  "
        f"max|d| {summary['max_abs']:.3e}  "
        f"up/down {summary['up']}/{summary['down']} z={summary['z']:+.1f}"
    )


def report_per_contig(runs: list[dict], baseline: dict, replicate: dict | None):
    print("=" * 74)
    print("3. PER-CONTIG PAIRED DELTAS  (vs baseline, matched on contig_name)")
    print("=" * 74)
    print("  Every contig is measured under both settings, so these are paired")
    print("  differences, not two independent samples. score is the product of")
    print("  the four components, so a component tells you which term moved.")
    print()

    names = sorted(set(baseline["contigs"]))
    keys = ("score",) + SCORE_COMPONENTS + ("coverage", "tpm")

    for run in runs:
        if run is baseline:
            continue
        tag = " (replicate = noise floor)" if run is replicate else ""
        print(f"  {run['label']}{tag}")
        if set(run["contigs"]) != set(baseline["contigs"]):
            print("    contig sets differ; skipped (see section 1)")
            print()
            continue
        for key in keys:
            if key not in next(iter(baseline["contigs"].values())):
                continue
            print_summary(
                key,
                paired_summary(column(baseline, names, key), column(run, names, key)),
            )
        print()

    if replicate is not None:
        print("  Read every row above against the replicate's row for the same")
        print("  metric. A setting whose 'changed' count and mean delta sit at or")
        print("  below the replicate's has no measurable effect on the score.")
    print()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "runs", nargs="+", metavar="LABEL=DIR",
        help="a transrate -o output directory, labelled",
    )
    parser.add_argument(
        "--baseline", default=None,
        help="label everything is compared against (default: the first run)",
    )
    parser.add_argument(
        "--replicate", default=None,
        help="label of a rerun at the baseline's settings; supplies the "
             "noise floor. Without it no difference can be called real",
    )
    args = parser.parse_args(argv)

    runs = []
    for item in args.runs:
        if "=" not in item:
            parser.error(f"expected LABEL=DIR, got {item!r}")
        label, _, directory = item.partition("=")
        runs.append(load_run(label, Path(directory)))

    by_label = {run["label"]: run for run in runs}
    if len(by_label) != len(runs):
        parser.error("labels must be unique")

    baseline = by_label[args.baseline] if args.baseline else runs[0]
    if args.baseline and args.baseline not in by_label:
        parser.error(f"--baseline {args.baseline} is not one of the runs")
    replicate = by_label.get(args.replicate) if args.replicate else None
    if args.replicate and replicate is None:
        parser.error(f"--replicate {args.replicate} is not one of the runs")

    print()
    print(f"assembly : {baseline['assembly'].get('assembly', '?')}")
    print(f"baseline : {baseline['label']}  ({baseline['dir']})")
    print(f"contigs  : {len(baseline['contigs'])}")
    print()

    valid = check_validity(runs, baseline)
    report_assembly_level(runs, baseline, replicate)
    report_per_contig(runs, baseline, replicate)
    return 0 if valid else 1


if __name__ == "__main__":
    sys.exit(main())
