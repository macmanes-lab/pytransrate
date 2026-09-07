#!/usr/bin/env python3
"""Which orthogroups change representative when the evaluator changes?

``compare_transrate_runs.py`` estimates how often two runs order a random
pair of contigs oppositely.  That is the right global measure, but ORP does
not compare random pairs: ``scripts/pick_best_contigs.py`` keeps the
highest-scoring member of each orthogroup, and orthogroup members are
near-duplicates whose scores sit close together.  This runs the real
selection against both runs' ``contigs.csv`` and reports exactly which
groups change winner.

Two things about where the inputs come from:

* The CSV must be the one ORP actually feeds the picker -- the
  **orthotransrate** run over ``merged.fasta`` (``oyster.py:580``), not a run
  over the finished ``.ORP.fasta``.  Scoring the final assembly answers a
  different question: those contigs have already survived selection.
* ``makeorthout`` deletes the ``*.groups`` files once it is done
  (``oyster.py:601``), so ``--orthogroups`` reads ``Orthogroups.txt``
  directly and rebuilds them in memory, exactly as ``makegroups`` does.

Usage::

    compare_orthogroup_picks.py old/contigs.csv new/contigs.csv \\
        --orthogroups .../Orthogroups.txt --out-prefix picks

    compare_orthogroup_picks.py old/contigs.csv new/contigs.csv \\
        --groups orthofuse_dir/
"""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import os
import sys
from collections import Counter

#: Column of contigs.csv holding the contig score. pick_best_contigs.py
#: hardcodes row[8]; output.py pins the same index. Kept as a constant so a
#: schema change breaks loudly here rather than silently picking a neighbour.
SCORE_COLUMN = 8

#: Groups whose two candidates differ by less than this are near-ties: the
#: two runs preferred different members of a pair they both consider
#: essentially equivalent. Well below the 6-decimal rounding of contigs.csv
#: would be meaningless, so this is a judgement about biology, not precision.
NEAR_TIE = 0.01


# ---------------------------------------------------------------------------
# SELECTION_RULE
#
# Mirrors pick_best_contigs.py exactly, and the details are load-bearing:
#
#   * ``max_score`` starts at 0.0 and the test is ``score > max_score``, so a
#     contig scoring exactly 0 is never selected and a group of them yields
#     no representative at all.
#   * Strictly greater, so ties keep the *first* member in file order -- the
#     rule the original awk one-liner had.
#   * A contig absent from the CSV is skipped rather than treated as zero.
#   * load_scores keeps the highest score when an id appears twice.
#
# Pass --pick-best to import the real implementation instead of this mirror,
# which removes any chance of the two drifting apart.
# ---------------------------------------------------------------------------


def load_scores(path: str) -> dict:
    if not os.path.isfile(path):
        sys.exit(f"contigs.csv not found at {path!r}")
    scores: dict = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)  # header
        for row in reader:
            if len(row) <= SCORE_COLUMN:
                continue
            try:
                score = float(row[SCORE_COLUMN])
            except ValueError:
                continue
            name = row[0]
            if name not in scores or score > scores[name]:
                scores[name] = score
    return scores


def best_in_group(members, scores):
    """The member ORP would keep, or None. See SELECTION_RULE."""
    best = 0.0
    want = None
    for name in members:
        score = scores.get(name)
        if score is not None and score > best:
            best = score
            want = name
    return want


def load_real_picker(path: str):
    """Import ORP's pick_best_contigs.py so the rule cannot drift."""
    spec = importlib.util.spec_from_file_location("pick_best_contigs", path)
    if spec is None or spec.loader is None:
        sys.exit(f"could not import {path!r}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.load_scores


def read_groups(groups_dir: str):
    """Yield ``(label, [members])`` from a directory of ``*.groups``."""
    paths = sorted(glob.glob(os.path.join(groups_dir, "*groups")))
    if not paths:
        sys.exit(f"no *.groups files under {groups_dir!r}")
    for path in paths:
        with open(path) as handle:
            members = [line.strip() for line in handle if line.strip()]
        yield os.path.basename(path), members


def read_orthogroups(path: str):
    """Yield ``(label, [members])`` from Orthogroups.txt.

    One group per line, ``OG0000001: contig_a contig_b ...``. The leading
    token is the group label and is dropped, as makegroups does.
    """
    with open(path) as handle:
        for index, line in enumerate(handle, start=1):
            tokens = line.split()
            if not tokens:
                continue
            yield f"{index}.groups", tokens[1:]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("old_csv", help="contigs.csv from the baseline run")
    parser.add_argument("new_csv", help="contigs.csv from the run being tested")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--groups", metavar="DIR", help="directory of *.groups")
    source.add_argument(
        "--orthogroups", metavar="FILE",
        help="Orthogroups.txt, when the *.groups files have been deleted",
    )
    parser.add_argument(
        "--out-prefix", metavar="P",
        help="write <P>.old.list, <P>.new.list and <P>.changed.tsv",
    )
    parser.add_argument(
        "--pick-best", metavar="PATH",
        help="ORP's scripts/pick_best_contigs.py, imported so the scoring "
             "rule is the real one rather than this script's mirror",
    )
    args = parser.parse_args(argv)

    reader = load_scores
    if args.pick_best:
        reader = load_real_picker(args.pick_best)
    old_scores = reader(args.old_csv)
    new_scores = reader(args.new_csv)

    groups = (
        read_groups(args.groups) if args.groups
        else read_orthogroups(args.orthogroups)
    )

    same = changed = only_old = only_new = neither = 0
    sizes: Counter = Counter()
    changed_sizes: Counter = Counter()
    near_ties = 0
    rows = []
    old_list, new_list = [], []

    for label, members in groups:
        if not members:
            continue
        size = len(members)
        sizes[size] += 1

        old_pick = best_in_group(members, old_scores)
        new_pick = best_in_group(members, new_scores)
        if old_pick:
            old_list.append(old_pick)
        if new_pick:
            new_list.append(new_pick)

        if old_pick is None and new_pick is None:
            neither += 1
        elif new_pick is None:
            only_old += 1
        elif old_pick is None:
            only_new += 1
        elif old_pick == new_pick:
            same += 1
        else:
            changed += 1
            changed_sizes[size] += 1
            # How decisively did each run prefer its own winner? A pair the
            # two runs both consider near-equivalent is a different finding
            # from one where each is confident and they disagree.
            old_margin = old_scores.get(old_pick, 0.0) - old_scores.get(new_pick, 0.0)
            new_margin = new_scores.get(new_pick, 0.0) - new_scores.get(old_pick, 0.0)
            if max(old_margin, new_margin) < NEAR_TIE:
                near_ties += 1
            rows.append((label, size, old_pick, new_pick, old_margin, new_margin))

    decided = same + changed
    total = decided + only_old + only_new + neither

    print()
    print("=" * 70)
    print("ORTHOGROUP REPRESENTATIVE CHANGES")
    print("=" * 70)
    print(f"  orthogroups                  {total:>10,}")
    print(f"  both runs picked a contig    {decided:>10,}")
    pct = (lambda n: f"   {n / decided:>7.2%}") if decided else (lambda n: "")
    print(f"    same representative        {same:>10,}{pct(same)}")
    print(f"    changed representative     {changed:>10,}{pct(changed)}")
    if changed:
        print(f"      of those, near-ties      {near_ties:>10,}"
              f"   {near_ties / changed:>7.2%}"
              f"  (both margins < {NEAR_TIE})")
    print(f"  only baseline picked         {only_old:>10,}")
    print(f"  only new picked              {only_new:>10,}")
    print(f"  neither picked               {neither:>10,}"
          "   (every member scored 0 or absent)")
    print()

    if changed:
        print("  By group size:")
        print(f"    {'members':>9}{'groups':>10}{'changed':>10}{'rate':>9}")
        for size in sorted(sizes):
            n, c = sizes[size], changed_sizes.get(size, 0)
            if n < 10:
                continue
            print(f"    {size:>9}{n:>10,}{c:>10,}{c / n:>9.2%}")
        print()
        print("  Largest disagreements (each run's margin over the other's pick):")
        rows.sort(key=lambda r: -max(r[4], r[5]))
        print(f"    {'group':<14}{'n':>4}  {'baseline margin':>16}{'new margin':>13}")
        for label, size, _o, _n, om, nm in rows[:10]:
            print(f"    {label:<14}{size:>4}  {om:>16.4f}{nm:>13.4f}")
        print()

    print(f"SUMMARY  groups={total}  decided={decided}  changed={changed}"
          f"  rate={changed / decided if decided else float('nan'):.4f}"
          f"  near_ties={near_ties}")

    if args.out_prefix:
        with open(f"{args.out_prefix}.old.list", "w") as handle:
            handle.write("".join(f"{name}\n" for name in old_list))
        with open(f"{args.out_prefix}.new.list", "w") as handle:
            handle.write("".join(f"{name}\n" for name in new_list))
        with open(f"{args.out_prefix}.changed.tsv", "w") as handle:
            handle.write("group\tmembers\tbaseline_pick\tnew_pick"
                         "\tbaseline_margin\tnew_margin\n")
            for label, size, old_pick, new_pick, om, nm in rows:
                handle.write(
                    f"{label}\t{size}\t{old_pick}\t{new_pick}\t{om:.6f}\t{nm:.6f}\n"
                )
        print(f"\nwrote {args.out_prefix}.old.list, .new.list, .changed.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
