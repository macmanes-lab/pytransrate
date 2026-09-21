#!/usr/bin/env python3
"""
build_manifest.py - scan the downloaded pair tree and write one TransRate job
per line.

Expects the layout fetch_pairs.py produces:

    $BASE/pairs/<pair_id>/assembly/*.fasta.gz
    $BASE/pairs/<pair_id>/reads/<RUN>_1.fastq.gz
    $BASE/pairs/<pair_id>/reads/<RUN>_2.fastq.gz

Writes a TSV: pair_id <tab> assembly <tab> left_csv <tab> right_csv
where left_csv/right_csv are comma-separated lists in matched order, which is
what TransRate expects for multi-run input.

Only pairs with exactly one R1/R2 set are emitted. Pairs with several runs are
skipped and listed, since combining them is a judgement call (concatenate? pick
one? assemble per tissue?) rather than something to guess at.

    python build_manifest.py --base /mnt/home/macmaneslab/macmanes/compare
    python build_manifest.py --base ... --skip-done          # omit finished pairs
    python build_manifest.py --base ... --allow-multi-run    # opt back in
"""

import argparse
import csv
import re
import sys
from pathlib import Path

ASM_EXT = (".fasta.gz", ".fasta", ".fa", ".fa.gz", ".fsa_nt.gz", ".fsa_nt", ".gz")


def find_assembly(pdir):
    adir = pdir / "assembly"
    if not adir.is_dir():
        return None
    cands = [p for p in sorted(adir.iterdir())
             if p.is_file() and p.name.endswith(ASM_EXT) and p.stat().st_size > 0]
    if not cands:
        return None
    # largest file wins if a directory somehow holds more than one
    return max(cands, key=lambda p: p.stat().st_size)


def find_pairs(pdir):
    """Return ([R1...], [R2...]) in matched order. Unpaired files are ignored."""
    rdir = pdir / "reads"
    if not rdir.is_dir():
        return [], [], []
    by_run = {}
    unpaired = []
    for p in sorted(rdir.iterdir()):
        if not p.is_file() or p.stat().st_size == 0:
            continue
        m = re.match(r"^(.+?)_([12])\.fastq(\.gz)?$", p.name)
        if m:
            by_run.setdefault(m.group(1), {})[m.group(2)] = p
        elif re.match(r"^(.+?)\.fastq(\.gz)?$", p.name):
            unpaired.append(p)          # ENA's orphan-read file; not usable here
    left, right = [], []
    for run in sorted(by_run):
        d = by_run[run]
        if "1" in d and "2" in d:
            left.append(d["1"])
            right.append(d["2"])
    return left, right, unpaired


def is_done(base, pair_id):
    """TransRate writes assemblies.csv into its output directory when it finishes."""
    return (base / "results" / pair_id / "assemblies.csv").exists()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path,
                    help="directory holding pairs/ (e.g. /mnt/home/macmaneslab/macmanes/compare)")
    ap.add_argument("--out", default="manifest.tsv", type=Path)
    ap.add_argument("--allow-multi-run", action="store_true",
                    help="also emit pairs with several runs, concatenated as "
                         "comma-separated lists. Off by default.")
    ap.add_argument("--skip-done", action="store_true",
                    help="omit pairs that already have results/<pair>/assemblies.csv")
    args = ap.parse_args()

    pairs_dir = args.base / "pairs"
    if not pairs_dir.is_dir():
        sys.exit(f"no pairs/ directory under {args.base}")

    jobs, skipped = [], []
    for pdir in sorted(p for p in pairs_dir.iterdir() if p.is_dir()):
        pid = pdir.name
        if args.skip_done and is_done(args.base, pid):
            skipped.append((pid, "already done"))
            continue
        asm = find_assembly(pdir)
        if not asm:
            skipped.append((pid, "no assembly file"))
            continue
        left, right, unpaired = find_pairs(pdir)
        if not left:
            why = ("only unpaired/single-end reads" if unpaired else "no read files")
            skipped.append((pid, why))
            continue
        if len(left) > 1 and not args.allow_multi_run:
            skipped.append((pid, f"{len(left)} read pairs (use --allow-multi-run)"))
            continue
        jobs.append((pid, str(asm),
                     ",".join(str(p) for p in left),
                     ",".join(str(p) for p in right)))

    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for j in jobs:
            w.writerow(j)

    print(f"{len(jobs)} jobs -> {args.out}")
    n_runs = sum(len(j[2].split(",")) for j in jobs)
    print(f"{n_runs} read pairs across all jobs")
    if skipped:
        print(f"\n{len(skipped)} pairs skipped:")
        for pid, why in skipped:
            print(f"  {pid}: {why}")
    if jobs:
        print(f"\nsubmit with:  sbatch --array=1-{len(jobs)}%20 run_transrate.sbatch")


if __name__ == "__main__":
    main()
