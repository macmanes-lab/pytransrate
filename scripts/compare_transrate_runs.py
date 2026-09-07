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

The same shape works for comparing *implementations* -- the Ruby transrate
ORP ships against this port, say.  Their CSVs carry identical columns in
identical order, and both round the same way, so nothing extra is needed.
Read it differently though: the port deliberately does not reproduce the
Ruby's scores (different aligner, different quantifier, a corrected soft-clip
cursor, an in-process replacement for salmon's postSample.bam), so the
question is not whether the numbers moved but whether the *ordering* did.
Section 5 answers that; ORP keeps the best-scoring member of each orthogroup,
so ordering is what actually reaches the assembly.

Every comparison also prints one grep-able SUMMARY line, so a sweep over
several datasets can be reduced to ``... | grep ^SUMMARY``.
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


#: Sequence-only values are compared numerically, not as text. Two
#: implementations of the same statistic can disagree in the last decimal --
#: the Ruby transrate computes prop_gc in a C extension, this port in Python
#: -- and a string comparison would call that a different assembly. Anything
#: above this is a real difference.
SEQUENCE_TOLERANCE = 1e-5


def values_differ(left: str, right: str) -> tuple[bool, float]:
    """Compare two CSV cells, numerically where both parse as numbers.

    Returns ``(differs, magnitude)``; magnitude is 0.0 for text comparisons.
    """
    if left == right:
        return False, 0.0
    a, b = as_float(left), as_float(right)
    if math.isnan(a) or math.isnan(b):
        return True, 0.0  # not numeric, and the text already differed
    delta = abs(a - b)
    return delta > SEQUENCE_TOLERANCE, delta


def check_validity(runs: list[dict], baseline: dict) -> bool:
    """The runs must have analysed the same assembly. Anything else invalidates
    the comparison, whether the runs differ by a mapping flag or by which
    implementation produced them."""
    print("=" * 74)
    print("1. VALIDITY  (sequence-only outputs must match)")
    print("=" * 74)

    ok = True
    base_names = set(baseline["contigs"])
    for run in runs:
        if run is baseline:
            continue
        problems = []
        rounding = []

        for key in ASSEMBLY_SEQUENCE_KEYS:
            if key == "assembly":
                continue  # the path differs by design
            if key not in baseline["assembly"]:
                continue
            differs, delta = values_differ(
                baseline["assembly"].get(key, ""), run["assembly"].get(key, "")
            )
            if differs:
                problems.append(
                    f"assemblies.csv {key}: "
                    f"{baseline['assembly'].get(key)} -> {run['assembly'].get(key)}"
                )
            elif delta:
                rounding.append(f"{key} ({delta:.1e})")

        names = set(run["contigs"])
        if names != base_names:
            problems.append(
                f"contig set differs: {len(base_names - names)} missing, "
                f"{len(names - base_names)} extra"
            )
        else:
            for key in CONTIG_SEQUENCE_KEYS:
                differing = 0
                worst = 0.0
                for name in base_names:
                    d, delta = values_differ(
                        baseline["contigs"][name][key], run["contigs"][name][key]
                    )
                    differing += d
                    worst = max(worst, delta)
                if differing:
                    problems.append(
                        f"contigs.csv {key}: {differing} contigs differ "
                        f"(max {worst:.3g})"
                    )
                elif worst:
                    rounding.append(f"{key} (max {worst:.1e})")

        if problems:
            ok = False
            print(f"  FAIL  {run['label']}")
            for problem in problems:
                print(f"          {problem}")
        else:
            print(f"  ok    {run['label']}")
            if rounding:
                print(
                    f"          within tolerance, last-decimal only: "
                    f"{', '.join(rounding)}"
                )

    if not ok:
        print(
            "\n  These runs did not analyse the same assembly, so the score\n"
            "  differences below are not what you think they are. Check that\n"
            "  every run was given the same FASTA."
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
        print("  No --replicate given, so there is no noise floor. Add one when")
        print("  comparing settings of one implementation, where the effect can")
        print("  be small enough for variance to matter. Comparing two different")
        print("  implementations it matters less -- those differences are")
        print("  structural, and section 5 is the one to read.")
    print()


# -- 3. score decomposition -------------------------------------------------


def contig_geomean(run: dict) -> float:
    """The geometric mean of the contig scores, recovered from the CSV.

    score.py computes the assembly score as

        score = geomean(contig scores) * (good_mappings / fragments)

    and assemblies.csv carries the score and both factors of the rate, so
    the geomean divides straight back out. Nothing is re-derived from
    contigs.csv: this is the number the run actually used.
    """
    rate = good_rate(run)
    score = as_float(run["assembly"].get("score"))
    if not rate or math.isnan(score):
        return math.nan
    return score / rate


def good_rate(run: dict) -> float:
    """``good_mappings / fragments``, the second factor of the score."""
    values = run["assembly"]
    rate = as_float(values.get("p_good_mapping"))
    if not math.isnan(rate) and rate:
        return rate
    fragments = as_float(values.get("fragments"))
    if not fragments or math.isnan(fragments):
        return math.nan
    return as_float(values.get("good_mappings")) / fragments


def report_decomposition(runs: list[dict], baseline: dict):
    """Split each score delta into its two exact causes.

    A score can rise because the contigs themselves scored better or because
    a larger share of fragments was called good -- different findings with
    different consequences, and the product form separates them exactly.
    """
    print("=" * 74)
    print("3. SCORE DECOMPOSITION")
    print("=" * 74)
    print("  score = geomean(contig scores) x (good_mappings / fragments)")
    print()

    if math.isnan(contig_geomean(baseline)):
        print("  No read metrics in this run, so there is no score to split.")
        print()
        return

    print(f"  {'run':<14}{'score':>12}{'contig geomean':>18}{'good/fragments':>18}")
    print("  " + "-" * 60)
    for run in runs:
        print(
            f"  {run['label']:<14}"
            f"{as_float(run['assembly']['score']):>12.5f}"
            f"{contig_geomean(run):>18.5f}"
            f"{good_rate(run):>18.5f}"
        )

    others = [run for run in runs if run is not baseline]
    if not others:
        print()
        return

    base_score = as_float(baseline["assembly"]["score"])
    base_geomean = contig_geomean(baseline)
    base_rate = good_rate(baseline)

    print()
    print("  Where each score delta comes from:")
    print()
    print(
        f"  {'run':<14}{'total':>11}{'from good-rate':>22}"
        f"{'from contig geomean':>24}"
    )
    print("  " + "-" * 70)
    for run in others:
        total = as_float(run["assembly"]["score"]) - base_score
        # score = G * r, so the delta splits as G_base*(r - r_base) for the
        # rate and (G - G_base)*r_base for the contigs, leaving only the
        # second-order cross term unaccounted for.
        from_rate = base_geomean * (good_rate(run) - base_rate)
        from_contigs = (contig_geomean(run) - base_geomean) * base_rate
        share = (
            lambda part: f" ({part / total:+.0%})" if abs(total) > 1e-12 else ""
        )
        print(
            f"  {run['label']:<14}{total:>+11.5f}"
            f"{f'{from_rate:+.5f}{share(from_rate)}':>22}"
            f"{f'{from_contigs:+.5f}{share(from_contigs)}':>24}"
        )
        residual = total - from_rate - from_contigs
        if abs(total) > 1e-12 and abs(residual) > 0.02 * abs(total):
            print(
                f"      unexplained {residual:+.5f} ({residual / total:+.0%}): "
                "CSV rounding plus the second-order cross term"
            )

    print()
    print("  A rise in the good-rate is more fragments being called good; a rise")
    print("  in the geomean is the contigs themselves scoring better. Read them")
    print("  against the accuracy components in the next section -- a score that")
    print("  climbs while sCnuc falls is not the assembly getting better.")
    print()


# -- 4. per contig ----------------------------------------------------------


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
    print("4. PER-CONTIG PAIRED DELTAS  (vs baseline, matched on contig_name)")
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


# -- 5. rank agreement ------------------------------------------------------

#: Random contig pairs drawn to estimate rank discordance. The estimate's
#: standard error is under 0.05 percentage points at this many draws, which
#: is far finer than any decision anyone makes from it.
RANK_PAIRS = 2_000_000


def rank_discordance(base: np.ndarray, other: np.ndarray, seed: int = 0) -> tuple:
    """How often the two runs disagree about which of two contigs is better.

    ORP's pick_best_contigs.py takes the highest-scoring member of each
    orthogroup, so what matters downstream is the *ordering* of contig
    scores, not their level. This samples random pairs and asks how often
    the two runs order them oppositely -- directly, how often a two-member
    orthogroup would change winner.

    Returns ``(discordance, tie_fraction)``. Pairs where either run scores
    the two contigs equally express no preference and are excluded from the
    discordance, since neither ordering is a disagreement.
    """
    n = base.size
    if n < 2:
        return math.nan, math.nan
    rng = np.random.default_rng(seed)  # fixed: the number must not drift
    i = rng.integers(0, n, size=RANK_PAIRS)
    j = rng.integers(0, n, size=RANK_PAIRS)
    keep = i != j
    i, j = i[keep], j[keep]

    delta_base = base[i] - base[j]
    delta_other = other[i] - other[j]
    decided = (delta_base != 0) & (delta_other != 0)
    n_decided = int(decided.sum())
    if not n_decided:
        return math.nan, 1.0
    discordant = int(
        np.sum(np.sign(delta_base[decided]) != np.sign(delta_other[decided]))
    )
    return discordant / n_decided, 1.0 - n_decided / delta_base.size


def classification(run: dict, names: list) -> np.ndarray | None:
    """Which contigs the run itself calls good, at its own optimal cutoff."""
    cutoff = as_float(run["assembly"].get("cutoff"))
    if math.isnan(cutoff):
        return None
    # Contig.classify keeps score >= cutoff; see CUTOFF_BOUNDARY in score.py.
    return column(run, names, "score") >= cutoff


def report_rank_agreement(runs: list[dict], baseline: dict):
    """Does the ordering survive, and would the same contigs be kept?"""
    print("=" * 74)
    print("5. RANK AGREEMENT  (would the same contigs be chosen?)")
    print("=" * 74)
    print("  ORP picks the best-scoring member of each orthogroup, so the")
    print("  ordering of contig scores is what reaches the assembly, not the")
    print("  scores themselves. Two runs can differ a lot in level and still")
    print("  make identical choices -- or agree closely and still not.")
    print()

    names = sorted(set(baseline["contigs"]))
    if "score" not in next(iter(baseline["contigs"].values())):
        print("  No contig scores in these runs (no reads).")
        print()
        return

    base_scores = column(baseline, names, "score")
    base_good = classification(baseline, names)

    for run in runs:
        if run is baseline:
            continue
        if set(run["contigs"]) != set(baseline["contigs"]):
            print(f"  {run['label']}: contig sets differ; skipped (see section 1)")
            continue

        other = column(run, names, "score")
        pearson = float(np.corrcoef(base_scores, other)[0, 1])
        disc, ties = rank_discordance(base_scores, other)

        print(f"  {run['label']}")
        print(f"    pearson r on contig score      {pearson:.5f}")
        print(
            f"    pairs ordered oppositely       {disc:.3%}"
            f"   (ties, no preference: {ties:.1%})"
        )

        good = classification(run, names)
        if base_good is not None and good is not None:
            both = int(np.sum(base_good & good))
            neither = int(np.sum(~base_good & ~good))
            lost = int(np.sum(base_good & ~good))
            gained = int(np.sum(~base_good & good))
            total = base_good.size
            print(
                f"    good/bad call at each run's own cutoff: "
                f"{(both + neither) / total:.2%} agree"
            )
            print(
                f"      good in both {both:>7}   bad in both {neither:>7}"
                f"   only baseline {lost:>6}   only {run['label']} {gained:>6}"
            )
        print()


def report_summary(runs: list[dict], baseline: dict):
    """One grep-able line per comparison, for looping over datasets."""
    names = sorted(set(baseline["contigs"]))
    has_scores = "score" in next(iter(baseline["contigs"].values()))
    base_scores = column(baseline, names, "score") if has_scores else None
    assembly = baseline["assembly"].get("assembly", "?")

    for run in runs:
        if run is baseline:
            continue
        fields = [
            f"{baseline['label']}->{run['label']}",
            f"d_score={as_float(run['assembly'].get('score')) - as_float(baseline['assembly'].get('score')):+.5f}",
        ]
        if not math.isnan(contig_geomean(baseline)):
            fields.append(
                f"d_goodrate={contig_geomean(baseline) * (good_rate(run) - good_rate(baseline)):+.5f}"
            )
            fields.append(
                f"d_geomean={(contig_geomean(run) - contig_geomean(baseline)) * good_rate(baseline):+.5f}"
            )
        if has_scores and set(run["contigs"]) == set(baseline["contigs"]):
            other = column(run, names, "score")
            disc, _ = rank_discordance(base_scores, other)
            fields.append(f"pearson={float(np.corrcoef(base_scores, other)[0, 1]):.5f}")
            fields.append(f"rank_disc={disc:.4f}")
        print(f"SUMMARY  {assembly}  " + "  ".join(fields))


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
    report_decomposition(runs, baseline)
    report_per_contig(runs, baseline, replicate)
    report_rank_agreement(runs, baseline)
    report_summary(runs, baseline)
    return 0 if valid else 1


if __name__ == "__main__":
    sys.exit(main())
