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

* The CSV must be the one ORP actually feeds the picker -- its
  **orthotransrate** step, scoring ``merged.fasta``, not the **transrate**
  step that scores the finished ``.ORP.fasta``.  Scoring the final assembly
  answers a different question: those contigs have already survived
  selection.
* ``--orthogroups`` reads ``Orthogroups.txt`` and rebuilds the groups in
  memory, which is what ORP itself now does.  Through ORP 3.x a
  ``makegroups`` step split that file into one ``<i>.groups`` file per
  orthogroup and ``makeorthout`` deleted them again on its way out; ORP
  4.0.0 removed that round-trip, and ``pick_best_contigs.py`` reads
  ``Orthogroups.txt`` directly, so the files are never written at all.
  ``--groups`` is therefore only useful against an archived pre-4.0.0 run
  directory.

Step names rather than ``oyster.py`` line numbers throughout: the line
numbers drifted within a single release.

Usage::

    compare_orthogroup_picks.py old/contigs.csv new/contigs.csv \\
        --orthogroups .../Orthogroups.txt --out-prefix picks

    compare_orthogroup_picks.py old/contigs.csv new/contigs.csv \\
        --groups archived_orthofuse_dir/      # pre-4.0.0 runs only
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

#: Below this share of group members present in a CSV, the comparison is
#: not measuring selection at all. See COVERAGE_CHECK.
MIN_MEMBER_COVERAGE = 0.90

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
# which removes any chance of the two drifting apart. From ORP 4.0.0 that
# covers the selection itself and not just the CSV loader -- see
# load_real_picker.
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
    """Import ORP's pick_best_contigs.py so the rule cannot drift.

    Returns ``(load_scores, best_in_group)``. Up to ORP 3.x the real
    ``best_in_group`` took a path to an ``<i>.groups`` file, so only the CSV
    loader could be reused and the selection rule stayed a mirror however
    this flag was passed. ORP 4.0.0 gave it the member list directly -- the
    signature this script's mirror already had -- so from that release the
    rule itself is importable. ``read_orthogroups`` exists only in the newer
    picker, which is what tells the two apart.
    """
    spec = importlib.util.spec_from_file_location("pick_best_contigs", path)
    if spec is None or spec.loader is None:
        sys.exit(f"could not import {path!r}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if hasattr(module, "read_orthogroups"):
        return module.load_scores, module.best_in_group
    return module.load_scores, best_in_group


def read_groups(groups_dir: str):
    """Yield ``(label, [members])`` from a directory of ``*.groups``.

    Only pre-4.0.0 ORP run directories have these; see the module docstring.
    """
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
    token is the group label and is dropped, as ORP does.

    Groups are labelled ``<i>.groups`` after the per-orthogroup files ORP
    wrote through 3.x, so a report reads the same whichever source it came
    from.
    """
    groups = []
    with open(path) as handle:
        for index, line in enumerate(handle, start=1):
            tokens = line.split()
            if not tokens:
                continue
            groups.append((f"{index}.groups", tokens[1:]))
    # Label order, not file order. ORP writes one line per group into
    # good.<run>.list in exactly this order -- lexicographic on the
    # <i>.groups name, so 1, 10, 100, 2, ... -- inherited from the glob its
    # picker used to do and kept deliberately, because that order reaches
    # cd-hit-est and so the final assembly. Matching it is what makes
    # --out-prefix's .old.list/.new.list diffable against a real
    # good.<run>.list, and makes a report from --orthogroups identical to one
    # from --groups, which sorts the same way by construction.
    groups.sort(key=lambda item: item[0])
    return groups


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("old_csv", help="contigs.csv from the baseline run")
    parser.add_argument("new_csv", help="contigs.csv from the run being tested")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--groups", metavar="DIR",
        help="directory of *.groups, from a pre-4.0.0 ORP run directory",
    )
    source.add_argument(
        "--orthogroups", metavar="FILE",
        help="Orthogroups.txt; the normal source, and the only one for a run "
             "made by ORP 4.0.0 or later",
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

    reader, picker = load_scores, best_in_group
    if args.pick_best:
        reader, picker = load_real_picker(args.pick_best)
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
    slots = found_old = found_new = 0

    for label, members in groups:
        if not members:
            continue
        size = len(members)
        sizes[size] += 1
        slots += size
        found_old += sum(1 for m in members if m in old_scores)
        found_new += sum(1 for m in members if m in new_scores)

        old_pick = picker(members, old_scores)
        new_pick = picker(members, new_scores)
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

    # COVERAGE_CHECK
    #
    # Both implementations floor a contig score at 0.01, so any contig
    # present in either CSV always beats best_in_group's 0.0 threshold. A
    # group nobody picks from is therefore a group whose members are simply
    # absent -- never one that scored zero. Absent members make groups
    # trivially agree, so low coverage reads as perfect reproduction unless
    # it is checked. The usual cause is scoring the finished assembly, which
    # holds one contig per group by construction: the one that already won.
    share_old = found_old / slots if slots else 0.0
    share_new = found_new / slots if slots else 0.0
    print(f"  group members                {slots:>10,}")
    print(f"    found in baseline CSV      {found_old:>10,}   {share_old:>7.2%}")
    print(f"    found in new CSV           {found_new:>10,}   {share_new:>7.2%}")
    if min(share_old, share_new) < MIN_MEMBER_COVERAGE:
        print()
        print("  *** most group members are missing from the scored assembly ***")
        print("  Groups whose members are absent cannot change representative, so")
        print("  the agreement below is measuring nothing. This is what scoring")
        print("  the finished assembly looks like -- it holds one contig per")
        print("  orthogroup already, the winner. Score the assembly the picker")
        print("  actually runs on (merged.fasta) with both implementations.")
    print()
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
          "   (no member present in either CSV)")
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
          f"  near_ties={near_ties}  member_coverage={min(share_old, share_new):.4f}")

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
