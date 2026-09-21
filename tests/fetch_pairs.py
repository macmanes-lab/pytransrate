#!/usr/bin/env python3
"""
fetch_pairs.py - download the assemblies and read sets listed in a pairs CSV,
keeping each assembly next to the reads it was built from.

On-disk layout (the point of this script):

    <dest>/
      manifest.csv              every file: pair_id, kind, path, bytes, md5, status
      pairs/
        tsa_GKPW/
          pair.json             the CSV row, so the directory is self-describing
          assembly/GKPW00000000.fsa_nt.gz
          reads/SRR26451323_1.fastq.gz
          reads/SRR26451323_2.fastq.gz
          ...

Given any file you can recover its pair from the path; given any pair_id you
can list both sides. manifest.csv is the flat view for downstream scripts.

Usage:
    python fetch_pairs.py --csv pairs.csv --dest /data/tx --dry-run
    python fetch_pairs.py --csv pairs.csv --dest /data/tx --jobs 6
    python fetch_pairs.py --csv pairs.csv --dest /data/tx --what assemblies
    python fetch_pairs.py --csv pairs.csv --dest /data/tx --only tsa_GKPW,tsa_GKPV
    python fetch_pairs.py --csv pairs.csv --dest /data/tx --max-runs 2
    python fetch_pairs.py --csv pairs.csv --dest /data/tx --verify

    pip install requests
"""

import argparse
import csv
import hashlib
import json
import re
import sys
import threading
import time
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FASTA_API = "https://www.ebi.ac.uk/ena/browser/api/fasta"
ENA_TSA_FTP = "https://ftp.ebi.ac.uk/pub/databases/ena/tsa/public"
ENA_WGS_FTP = "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public"
NCBI_TSA_FTP = "https://ftp.ncbi.nlm.nih.gov/sra/wgs_aux"
NCBI_WGS = "https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux"
CHUNK = 1 << 20
UA = {"User-Agent": "transcriptome-pair-fetcher/1.0"}

_manifest_lock = threading.Lock()
_print_lock = threading.Lock()

# running totals, for the progress line
_stats = {"files_done": 0, "files_total": 0, "bytes": 0,
          "pairs_done": 0, "pairs_total": 0, "start": None}


def log(msg, level="  "):
    """Timestamped, thread-safe line."""
    with _print_lock:
        print(f"[{datetime.now():%H:%M:%S}] {level}{msg}", flush=True)


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:,.1f} {unit}" if unit != "B" else f"{n:,.0f} B"
        n /= 1024


def hms(seconds):
    return str(timedelta(seconds=int(seconds)))


def progress():
    """One-line running total, printed after each file."""
    el = time.time() - _stats["start"]
    rate = _stats["bytes"] / el if el > 0 else 0
    done, total = _stats["files_done"], _stats["files_total"]
    eta = ""
    if done and total and done < total:
        eta = f" | eta {hms((total - done) * (el / done))}"
    log(f"progress: {done}/{total} files, "
        f"{_stats['pairs_done']}/{_stats['pairs_total']} pairs, "
        f"{human(_stats['bytes'])}, {hms(el)} elapsed, "
        f"{human(rate)}/s{eta}", level="   ")


def rm(path):
    """Path.unlink(missing_ok=True) needs Python 3.8; this works on 3.6+."""
    try:
        path.unlink()
    except OSError:
        pass


def md5sum(path):
    h = hashlib.md5()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(CHUNK), b""):
            h.update(b)
    return h.hexdigest()


def gz_ok(path):
    """Cheap integrity check: gzip magic + a full decompress pass."""
    import gzip
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(CHUNK):
                pass
        return True
    except Exception:
        return False


def download(url, dest, session, expect_md5=None):
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size > 0:
        if expect_md5:
            if md5sum(dest) == expect_md5:
                return "cached", dest.stat().st_size
            dest.unlink()
        elif dest.suffix == ".gz":
            if gz_ok(dest):
                return "cached", dest.stat().st_size
            dest.unlink()
        else:
            return "cached", dest.stat().st_size

    headers = dict(UA)
    mode = "wb"
    if dest.exists():
        headers["Range"] = f"bytes={dest.stat().st_size}-"
        mode = "ab"

    with session.get(url, headers=headers, stream=True, timeout=300) as r:
        if r.status_code == 416:
            return "cached", dest.stat().st_size
        r.raise_for_status()
        if r.status_code == 200:
            mode = "wb"
        with dest.open(mode) as fh:
            for chunk in r.iter_content(CHUNK):
                fh.write(chunk)

    size = dest.stat().st_size
    if expect_md5 and md5sum(dest) != expect_md5:
        return "md5-mismatch", size
    return "downloaded", size


def ena_fastqs(run, session, platform="ILLUMINA"):
    """[(url, filename, md5), ...] for a run accession. Handles SE and PE.

    A study or BioProject accession expands to all of its runs, so the platform
    gate matters here too: it keeps 454/PacBio/Nanopore runs out of a mixed
    project rather than downloading them and finding out later.
    """
    params = {"accession": run, "result": "read_run",
              "fields": ("run_accession,fastq_ftp,fastq_md5,"
                         "instrument_platform,instrument_model"),
              "format": "tsv"}
    r = session.get(ENA_FILEREPORT, params=params, headers=UA, timeout=90)
    r.raise_for_status()
    lines = [l for l in r.text.splitlines() if l.strip()]
    if len(lines) < 2:
        return []
    hdr = lines[0].split("\t")
    out = []
    for line in lines[1:]:
        d = dict(zip(hdr, line.split("\t")))
        plat = (d.get("instrument_platform") or "").upper()
        if platform and plat and plat != platform.upper():
            continue
        paths = (d.get("fastq_ftp") or "").split(";")
        md5s = (d.get("fastq_md5") or "").split(";")
        for i, p in enumerate(paths):
            if not p:
                continue
            out.append((f"https://{p}", Path(p).name,
                        md5s[i] if i < len(md5s) and md5s[i] else None))
    return out


def assembly_urls(accession, provider="ncbi"):
    """Candidate download URLs for a TSA/WGS set accession like GAAA00000000.

    IMPORTANT: the 12-character master accession (GAAA01000000) identifies the
    record that *describes* the set, not the sequences. Requesting FASTA for it
    returns a near-empty file (~20 bytes). The sequences live under the
    6-character set prefix (GAAA01), so that is what we ask for:

      ENA API      .../ena/browser/api/fasta/GAAA01?download=true&gzip=true
      ENA FTP TSA  .../ena/tsa/public/gaa/GAAA01.fasta.gz
      ENA FTP WGS  .../ena/wgs/public/gaa/GAAA01.fasta.gz
      NCBI FTP     .../sra/wgs_aux/GA/AA/GAAA01/GAAA01.1.fsa_nt.gz

    TSA and WGS have separate FTP trees, and a transcriptome can appear under
    either, so both are tried. Every candidate is validated after download (see
    validate_assembly) and a file that fails validation moves on to the next.
    """
    m = re.match(r"^([A-Z]{4,6})(\d{2})?", accession.upper())
    if not m:
        return []
    prefix = m.group(1)
    ver = m.group(2)
    versions = [ver] if ver and ver != "00" else []
    versions += [v for v in ("01", "02", "03") if v not in versions]

    ena, ncbi = [], []
    for v in versions:
        stem = f"{prefix}{v}"                 # GBIS01 - the set
        a, b = prefix[:2], prefix[2:4]
        # NCBI serves TSA sets from a stable, documented layout. These
        # accessions are GenBank TSA, so this is the most reliable route.
        ncbi.append(f"{NCBI_TSA_FTP}/{a}/{b}/{stem}/{stem}.1.fsa_nt.gz")
        ncbi.append(f"{NCBI_WGS}/{a}/{b}/{stem}/{stem}.1.fsa_nt.gz")
        # ENA keeps set files on FTP, under separate tsa/ and wgs/ trees. The
        # subdirectory is the lowercased prefix; 3 and 2 characters are both
        # seen, so try both.
        for sub in (prefix[:3].lower(), prefix[:2].lower()):
            ena.append(f"{ENA_TSA_FTP}/{sub}/{stem}.fasta.gz")
            ena.append(f"{ENA_WGS_FTP}/{sub}/{stem}.fasta.gz")
        # Browser API last: it returns the master record (~20 B) for set
        # accessions, so it is kept only as a long shot, never as the default.
        ena.append(f"{ENA_FASTA_API}/{stem}?download=true&gzip=true")
    return ena + ncbi if provider == "ena" else ncbi + ena


def validate_assembly(path, min_bytes=2048):
    """(ok, reason). Catches the master-record-only file and other empties.

    Reads the head of the file rather than all of it, so this stays cheap on a
    multi-GB assembly.
    """
    import gzip
    try:
        size = path.stat().st_size
    except OSError:
        return False, "missing"
    if size < min_bytes:
        return False, f"too small ({size} B) - probably the master record, not the set"

    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            head = fh.read(200_000)
    except Exception as e:
        return False, f"unreadable ({type(e).__name__})"

    if ">" not in head:
        return False, "no FASTA header found"
    seq_chars = sum(len(l.strip()) for l in head.splitlines()
                    if l and not l.startswith(">"))
    if seq_chars < 100:
        return False, f"header present but only {seq_chars} sequence characters"
    return True, f"{head.count('>')}+ records in first {len(head)//1024} KB"


def handle_pair(row, dest, args, session, manifest):
    pid = row["pair_id"]
    pdir = dest / "pairs" / pid
    pdir.mkdir(parents=True, exist_ok=True)
    (pdir / "pair.json").write_text(json.dumps(row, indent=2), encoding="utf-8")

    jobs = []   # (kind, url, path, md5)

    # ---- assembly
    if args.what in ("all", "assemblies"):
        url = (row.get("assembly_url") or "").strip()
        acc = (row.get("assembly_accession") or "").strip()
        if url:
            name = url.split("?")[0].rstrip("/").split("/")[-1] or f"{acc}.fasta.gz"
            if not name.endswith((".gz", ".fasta", ".fa", ".fsa_nt")):
                name = f"{acc or name}.fasta.gz"
            jobs.append(("assembly", url, pdir / "assembly" / name, None))
        elif row.get("assembly_source") == "tsa" and acc:
            cands = assembly_urls(acc, args.provider)
            if cands:
                jobs.append(("assembly", cands,
                                     pdir / "assembly" / f"{acc}.fasta.gz", None))
        else:
            with _manifest_lock:
                manifest.append({"pair_id": pid, "kind": "assembly", "path": "",
                                 "bytes": 0, "md5": "",
                                 "status": "no-url (manual: see notes column)"})

    # ---- reads
    if args.what in ("all", "reads"):
        runs = [r for r in (row.get("run_accessions") or "").split(";") if r]
        if args.max_runs:
            runs = runs[:args.max_runs]
        for run in runs:
            if args.dry_run:
                jobs.append(("reads", f"<ENA lookup: {run}>",
                             pdir / "reads" / f"{run}_[12].fastq.gz", None))
                continue
            try:
                files = ena_fastqs(run, session, args.platform)
            except Exception as e:
                with _manifest_lock:
                    manifest.append({"pair_id": pid, "kind": "reads", "path": run,
                                     "bytes": 0, "md5": "",
                                     "status": f"ena-lookup-failed: {e}"})
                continue
            if not files:
                with _manifest_lock:
                    manifest.append({"pair_id": pid, "kind": "reads", "path": run,
                                     "bytes": 0, "md5": "",
                                     "status": f"no {args.platform or 'RNA-seq'} "
                                               "fastq in ENA (try prefetch)"})
            for url, name, md5 in files:
                jobs.append(("reads", url, pdir / "reads" / name, md5))

    if not args.dry_run:
        with _manifest_lock:
            _stats["files_total"] += len(jobs)
        log(f"[{pid}] {len(jobs)} files queued "
            f"({sum(1 for j in jobs if j[0] == 'assembly')} assembly, "
            f"{sum(1 for j in jobs if j[0] == 'reads')} reads) - {row.get('organism', '')}")

    results = []
    for kind, url, path, md5 in jobs:
        if args.dry_run:
            u = url[0] if isinstance(url, list) else url
            log(f"[{pid}] {kind:9} {u}\n              -> {path}")
            continue

        urls = url if isinstance(url, list) else [url]
        log(f"[{pid}] START    {kind:9} {path.name}")
        t0 = time.time()
        status, size, last = "failed", 0, ""
        for u in urls:
            try:
                status, size = download(u, path, session, md5)
            except Exception as e:
                last = f"{type(e).__name__}: {e}"
                continue
            # an assembly that downloads but has no sequence in it is a failure,
            # not a success - fall through to the next candidate URL
            if kind == "assembly":
                ok, why = validate_assembly(path, args.min_assembly_bytes)
                if not ok:
                    log(f"[{pid}] REJECT   {kind:9} {path.name}  {why}\n"
                        f"              url: {u}", level="!! ")
                    last = why
                    rm(path)
                    status, size = "failed", 0
                    continue
                log(f"[{pid}] VALID    {kind:9} {path.name}  {why}")
            break
        dt = time.time() - t0

        if status == "failed":
            status = f"failed ({last})"
            log(f"[{pid}] FAILED   {kind:9} {path.name}  {last}", level="!! ")
        elif status == "cached":
            log(f"[{pid}] CACHED   {kind:9} {path.name}  {human(size)}")
        elif status.startswith("MD5"):
            log(f"[{pid}] BAD MD5  {kind:9} {path.name}", level="!! ")
        else:
            rate = size / dt if dt > 0 else 0
            log(f"[{pid}] DONE     {kind:9} {path.name}  "
                f"{human(size)} in {hms(dt)} ({human(rate)}/s)")

        results.append((kind, path, size, md5, status))
        with _manifest_lock:
            manifest.append({"pair_id": pid, "kind": kind, "path": str(path),
                             "bytes": size, "md5": md5 or "", "status": status})
            _stats["files_done"] += 1
            _stats["bytes"] += size
        progress()

    return pid, results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="pairs.csv")
    ap.add_argument("--dest", default="/data/transcriptomes", type=Path)
    ap.add_argument("--what", choices=["all", "assemblies", "reads"], default="all")
    ap.add_argument("--only", help="comma-separated pair_ids")
    ap.add_argument("--max-runs", type=int, default=0,
                    help="cap runs per pair (0 = all). Useful for a sampling pass.")
    ap.add_argument("--provider", choices=["ncbi", "ena"], default="ncbi",
                    help="which archive to try first for ASSEMBLIES. Default ncbi: "
                         "ENA's browser API returns the master record rather than "
                         "the sequence set for TSA accessions, so NCBI's wgs_aux "
                         "layout is the dependable route. Reads always come from "
                         "ENA regardless.")
    ap.add_argument("--min-assembly-bytes", type=int, default=2048,
                    help="reject an assembly smaller than this and try the next "
                         "source. Guards against ENA returning the master record "
                         "instead of the sequence set.")
    ap.add_argument("--platform", default="ILLUMINA",
                    help="only download runs on this ENA instrument_platform. "
                         "Default ILLUMINA. Pass --platform '' to allow any.")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verify", action="store_true",
                    help="re-check existing files against manifest md5s, download nothing")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        sys.exit(f"CSV not found: {csv_path}")
    with csv_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    if args.only:
        keep = {x.strip() for x in args.only.split(",")}
        rows = [r for r in rows if r["pair_id"] in keep]
    if not rows:
        sys.exit("no pairs selected")

    args.dest.mkdir(parents=True, exist_ok=True)
    manifest_path = args.dest / "manifest.csv"

    if args.verify:
        if not manifest_path.exists():
            sys.exit("no manifest.csv to verify against")
        bad = 0
        with manifest_path.open(newline="", encoding="utf-8") as fh:
            for m in csv.DictReader(fh):
                p = Path(m["path"])
                if not m["path"] or not p.exists():
                    print(f"MISSING  {m['pair_id']}  {m['path']}")
                    bad += 1
                elif m["md5"] and md5sum(p) != m["md5"]:
                    print(f"BAD MD5  {m['pair_id']}  {p.name}")
                    bad += 1
                elif m["kind"] == "assembly":
                    ok, why = validate_assembly(p)
                    if not ok:
                        print(f"BAD ASM  {m['pair_id']}  {p.name}  {why}")
                        bad += 1
        print(f"\n{'all files verified' if not bad else f'{bad} problems'}")
        return

    _stats["start"] = time.time()
    _stats["pairs_total"] = len(rows)
    started = datetime.now()

    print("=" * 72)
    print(f"fetch_pairs.py  started {started:%Y-%m-%d %H:%M:%S}")
    print(f"  source csv : {csv_path}")
    print(f"  destination: {args.dest}")
    print(f"  pairs      : {len(rows)}")
    print(f"  fetching   : {args.what}")
    print(f"  platform   : {args.platform or 'any'}")
    print(f"  parallel   : {args.jobs}")
    if args.max_runs:
        print(f"  max runs   : {args.max_runs} per pair")
    print("=" * 72, flush=True)

    session = requests.Session()
    manifest = []
    failures = []

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(handle_pair, r, args.dest, args, session, manifest): r
                for r in rows}
        for i, fut in enumerate(as_completed(futs), 1):
            row = futs[fut]
            try:
                pid, results = fut.result()
                ok = sum(1 for r in results if r[4] in ("downloaded", "cached"))
                tot = sum(r[2] for r in results)
                with _manifest_lock:
                    _stats["pairs_done"] += 1
                bad = [r for r in results if r[4] not in ("downloaded", "cached")]
                log(f"[{pid}] PAIR COMPLETE  {ok}/{len(results)} files, {human(tot)}"
                    + (f", {len(bad)} problem(s)" if bad else ""),
                    level="** ")
                for r in bad:
                    failures.append((pid, r[1].name, r[4]))
            except Exception as e:
                log(f"[{row['pair_id']}] PAIR FAILED  {e}", level="!! ")
                failures.append((row["pair_id"], "-", str(e)))

    if args.dry_run:
        return

    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair_id", "kind", "path", "bytes",
                                           "md5", "status"])
        w.writeheader()
        w.writerows(manifest)

    total = sum(m["bytes"] for m in manifest)
    elapsed = time.time() - _stats["start"]
    n_asm = sum(1 for m in manifest if m["kind"] == "assembly"
                and m["status"] in ("downloaded", "cached"))
    n_rd = sum(1 for m in manifest if m["kind"] == "reads"
               and m["status"] in ("downloaded", "cached"))

    print("\n" + "=" * 72)
    print("DONE")
    print(f"  started    : {started:%Y-%m-%d %H:%M:%S}")
    print(f"  finished   : {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"  elapsed    : {hms(elapsed)}")
    print(f"  pairs      : {_stats['pairs_done']}/{len(rows)}")
    print(f"  assemblies : {n_asm}")
    print(f"  read files : {n_rd}")
    print(f"  downloaded : {human(total)}"
          f"  (avg {human(total / elapsed if elapsed else 0)}/s)")
    print(f"  manifest   : {manifest_path}")
    print("=" * 72)

    if failures:
        print(f"\n{len(failures)} problem(s):")
        for pid, name, why in failures[:30]:
            print(f"  {pid:<20} {name:<34} {why}")
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more (see {manifest_path})")
        print("\nRe-run the same command to retry; completed files are skipped.")
        sys.exit(1)
    print("\nAll files retrieved. Next: ./submit_all.sh")


if __name__ == "__main__":
    main()
