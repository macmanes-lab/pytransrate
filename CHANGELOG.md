# Changelog

Notable changes to pytransrate. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

Versions up to and including 1.0.3 are the original Ruby
[transrate](https://github.com/blahah/transrate). 2.0.0 is a rewrite in
Python and **does not reproduce their scores** — see *Changed* below.

## [2.0.0] — 2026-09-07

First release of the Python implementation. The Ruby version is pinned to
Ruby 2.2.0 (end of life since 2018), ships as an x86_64-linux-only bundle,
and depends on three abandoned binaries; none of that could be carried
forward, and updating the dependencies changed real behaviour.

### Changed — scores are not comparable to 1.0.3

The four score components, the contig score and the assembly score are
computed by the same formulas, but their inputs changed:

- **snap-aligner 1.0dev.96 → 2.0.5.** SNAP 2.0 introduced soft clipping, on
  by default. The old `bam-read` helper advanced its reference cursor on
  soft clips as though they consumed reference, which they do not, shifting
  a clipped read's whole coverage contribution rightward. That was inert
  under an aligner that never clipped and is corrupting under one that does.
  Coverage now follows the SAM spec and agrees with `samtools depth`
  (`SOFT_CLIP_FIX`).

  Measured against 1.0.3 on three assemblies of one library: `bases_uncovered`
  rises 3.1–4.5×, `sCnuc` rises, `sCcov`/`sCord`/`sCseg` fall, and the
  assembly score drops 0.008–0.070. The good-mapping *rate* is unchanged to
  ±0.001, so the two implementations agree about the reads and disagree about
  the contigs. 14.5–19% of contig pairs are ordered oppositely, so
  score-ranked downstream selection will not make identical choices.

- **salmon 0.8.2 → 2.7.0.** `--useErrorModel` is now a hard error and
  `--sampleOut` is accepted but does nothing, so the `postSample.bam` the old
  pipeline consumed is never written.

- **Fragment assignment is a redesign, not a port** (`ASSIGNMENT_MODEL`).
  salmon 0.8.2 sampled from its own posterior, which is not recoverable from
  outside. Assignment is now a deterministic maximum-a-posteriori choice
  using salmon's abundances as the prior and alignment edit distance as the
  likelihood, streamed straight into metric accumulation rather than staged
  through a BAM. Ties break on prior then name, never on BAM order.

- **`bam-read` is gone.** Per-contig accumulation is in-process.

- **The default snap flags differ from the Ruby's**
  (`MULTI_ALIGNMENT_SETTINGS`). `-H 300000 -D 5 -om 5 -omax 10` makes snap
  2.0.5 die with SIGFPE on real assemblies. Defaults are now `-H 4000 -D 2
  -om 2 -omax 10 -mpc 1`.

- **`-mcp` is no longer passed** (`MAX_CANDIDATE_POOL`). The Ruby's
  `-mcp 10000000000000` overflowed snap's `atoi()` to an arbitrary value —
  measured at 1316134912 — so it never meant what it looked like. Pass
  `--max-candidate-pool` explicitly to restore it.

- **The run report goes to stdout**, the banner to stderr.

Preserved deliberately, because there is no oracle to call a change an
improvement: `BINNING_QUIRK`, `FRAGMENT_ESTIMATOR`, `CUTOFF_BOUNDARY`,
`STATS_QUIRK`, `ORF_CASE_SENSITIVITY`.

The CSV column order is unchanged from 1.0.3 in both files. Downstream tools
read them positionally and `tests/test_output.py` pins the indices.

### Added

- `--max-alignments-per-contig`, `--max-alignments-per-pair`,
  `--multi-edit-distance`, `--extra-search-depth`, `--max-seed-hits`,
  `--edit-distance`, `--max-candidate-pool`, `--no-error-model`,
  `--location-size`, `--seed-size`, `--keep-bam`, `--no-banner`,
  `--loglevel`, `--version`.
- The invocation is logged as the first line of every run, shell-quoted, so a
  log identifies the settings that produced it.
- Contig extremes and the complete read-metrics block are reported during the
  run, formatted exactly as `assemblies.csv` rounds them.
- Soft clipping is counted and reported per run — clipped alignments, clipped
  bases, leading clipped bases — since the BAM it comes from is deleted by
  default.
- `scripts/compare_transrate_runs.py`: compares runs from their CSVs, with a
  sequence-only validity control, an exact split of any score difference into
  its good-rate and contig-geomean causes, per-contig paired deltas, and rank
  agreement between two runs.
- `scripts/compare_orthogroup_picks.py`: runs the Oyster River Protocol's
  orthogroup selection against two runs and reports exactly which groups
  change representative.

### Fixed

- Records that snap-aligner writes with a CIGAR disagreeing with the read
  length — a long-standing snap bug at contig boundaries, present in every
  2.0.x release — killed a run mid-scoring with a bare `OSError`. They are
  now skipped and counted, because htslib consumes the record before it
  validates and the stream is already positioned at the next one. A run of
  100 consecutive failures is still fatal: that is a truncated BAM, where the
  stream really is desynchronised (`MALFORMED_RECORDS`).
- Coverage no longer advances the reference cursor on soft clips
  (`SOFT_CLIP_FIX`).
- A failed snap index build no longer leaves a partial index that a later run
  would trust; the marker file is checked, not the directory.
- snap exiting 0 after failing to write a usable BAM is now caught.

### Performance

Profiled, not guessed. On 26,950 contigs and 21,390 fragments, 21.0s → 12.2s
end to end, with byte-identical scores:

| | before | after | |
| --- | --- | --- | --- |
| `_log_segment_likelihood` | 4.57s | 0.10s | 44.6× |
| `bin_coverage` | 3.00s | 0.37s | 8.1× |
| `prob_not_segmented` | 6.03s | 2.27s | 2.7× |
| `longest_orf` | 4.84s | 2.69s | 1.8× |

### Removed

- `--install-deps`. It drove the `bindeps` gem, fetching binaries from URLs
  that no longer resolve. snap-aligner and salmon are ordinary bioconda
  packages.
- `--reference` and all reference-based (CRB-BLAST) metrics. The Ruby routed
  these through the unmaintained `crb-blast` gem. The flag is still parsed
  and raises, so a script using it gets an error rather than silently
  different output.

[2.0.0]: https://github.com/macmanes-lab/transrate/releases/tag/v2.0.0
