#!/usr/bin/env python3
"""
collect_results.py - merge every finished TransRate run into one CSV carrying
all of its assembly metrics.

TransRate writes <output>/assemblies.csv with one row per assembly and a wide
set of metric columns (n_seqs, smallest, largest, n_bases, mean_orf_percent,
n50, score, optimal_score, cutoff, weighted, p_good, p_bases_covered,
fragments_mapped, good_mappings, bad_mappings, and so on). Column sets can
differ slightly between runs, so this takes the union and leaves blanks rather
than dropping columns or misaligning them.

    python collect_results.py --base /mnt/home/macmaneslab/macmanes/compare
    python collect_results.py --base ... --out metrics.csv
    python collect_results.py --base ... --contigs        # also per-contig files

Output columns: pair_id, organism, lineage, platform, assembly_accession,
n_read_pairs, then every TransRate metric found.
"""

import argparse
import csv
import sys
from pathlib import Path

# stable leading columns; everything else is appended in first-seen order
LEAD = ["pair_id", "organism", "lineage", "platform", "assembly_accession",
        "n_read_pairs", "transrate_dir"]


def load_pair_metadata(path):
    """pair_id -> row from pairs_seed.csv / pairs.csv, if the file is there."""
    if not path or not Path(path).exists():
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["pair_id"]: r for r in csv.DictReader(fh)}


def load_manifest(path):
    """pair_id -> number of read pairs, from the job manifest."""
    out = {}
    if not path or not Path(path).exists():
        return out
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3:
                out[row[0]] = len(row[2].split(","))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/mnt/home/macmaneslab/macmanes/compare",
                    type=Path)
    ap.add_argument("--out", type=Path,
                    help="default: <base>/transrate_metrics.csv")
    ap.add_argument("--pairs-csv", type=Path,
                    help="default: <base>/pairs_seed.csv, for organism/platform columns")
    ap.add_argument("--contigs", action="store_true",
                    help="also merge per-contig files into transrate_contigs.csv "
                         "(large: one row per contig across all assemblies)")
    args = ap.parse_args()

    base = args.base
    out = args.out or base / "transrate_metrics.csv"
    meta = load_pair_metadata(args.pairs_csv or base / "pairs_seed.csv")
    nreads = load_manifest(base / "manifest.tsv")

    results = base / "results"
    if not results.is_dir():
        sys.exit(f"no results/ directory under {base}")

    rows, metric_cols = [], []
    for adir in sorted(p for p in results.iterdir() if p.is_dir()):
        pid = adir.name
        acsv = adir / "assemblies.csv"
        if not acsv.exists() or acsv.stat().st_size == 0:
            continue
        with acsv.open(newline="", encoding="utf-8") as fh:
            for rec in csv.DictReader(fh):
                m = meta.get(pid, {})
                row = {
                    "pair_id": pid,
                    "organism": m.get("organism", ""),
                    "lineage": m.get("lineage", ""),
                    "platform": m.get("platform", ""),
                    "assembly_accession": m.get("assembly_accession", ""),
                    "n_read_pairs": nreads.get(pid, ""),
                    "transrate_dir": str(adir),
                }
                for k, v in rec.items():
                    if k is None:
                        continue
                    key = k.strip()
                    if key and key not in metric_cols and key not in LEAD:
                        metric_cols.append(key)
                    row[key] = v
                rows.append(row)

    if not rows:
        sys.exit(f"no assemblies.csv found under {results} - nothing finished yet")

    header = LEAD + metric_cols
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in header})

    print(f"{len(rows)} assemblies, {len(metric_cols)} metrics -> {out}")
    print("metrics: " + ", ".join(metric_cols))

    # which jobs have not landed
    manifest = base / "manifest.tsv"
    if manifest.exists():
        want = {r[0] for r in csv.reader(manifest.open(newline=""), delimiter="\t") if r}
        got = {r["pair_id"] for r in rows}
        missing = sorted(want - got)
        if missing:
            print(f"\n{len(missing)} of {len(want)} jobs have no results yet:")
            for p in missing:
                print(f"  {p}")

    if args.contigs:
        cout = base / "transrate_contigs.csv"
        n, hdr = 0, None
        with cout.open("w", newline="", encoding="utf-8") as fh:
            w = None
            for adir in sorted(p for p in results.iterdir() if p.is_dir()):
                for ccsv in sorted(adir.rglob("contigs.csv")):
                    with ccsv.open(newline="", encoding="utf-8") as cf:
                        rd = csv.DictReader(cf)
                        for rec in rd:
                            if w is None:
                                hdr = ["pair_id"] + [c for c in rd.fieldnames if c]
                                w = csv.DictWriter(fh, fieldnames=hdr,
                                                   extrasaction="ignore")
                                w.writeheader()
                            rec["pair_id"] = adir.name
                            w.writerow(rec)
                            n += 1
        print(f"\n{n} contigs -> {cout}")


if __name__ == "__main__":
    main()
