# Changelog

Notable changes to pytransrate. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

Versions up to and including 1.0.3 are the original Ruby
[transrate](https://github.com/blahah/transrate). 2.0.0 is a rewrite in
Python and **does not reproduce their scores** — see *Changed* below.

## [2.1.0] — 2026-09-09

A performance release. `-t/--threads` now applies to scoring as well as to
mapping and quantification, and the serial work each process does is 1.9×
faster per alignment and 5.3× faster per contig. **Scores are unchanged**
apart from `p_seq_true`, which moved by up to 3.3e-15 when it was made exact
so that a run no longer depends on its thread count — see *Changed*.

### Added

- `-t/--threads` now divides the scoring step as well as snap and salmon.
  Fragments are independent units of work and every accumulator they feed is
  additive, so the BAM is split by striding: each worker reads the whole file
  and takes the fragments where `index % workers == worker`, writing into
  anonymous mmaps shared through `fork` (`STRIDING`). Nothing large is
  pickled and nothing comes back through a pipe.

  The speedup is bounded by the share each worker pays regardless —
  decompression and a fragment-boundary check on every record, ~0.26µs of
  the ~1.5µs a record costs: ~3.7× at 8 workers, ~5.2× at 32, asymptotically
  ~6×. The alternative of scanning for byte offsets first was measured and
  comes out level, while needing an index the BAM does not have.

  Falls back to one process wherever `fork` is unavailable, rather than
  pretending. A worker killed by the OOM killer is now reported as such,
  with the per-worker footprint, instead of hanging the parent on a queue.
- Coverage integration and segmentation scoring are divided the same way,
  by contiguous ranges of contigs, on assemblies past 5,000 of them.

### Changed

- **`p_seq_true` is exact, and no longer depends on `--threads`.** It was the
  one accumulator carried as a running float, and float addition is not
  associative, so dividing a contig's reads between workers moved it by up to
  3.3e-15 — enough to reach the sixth decimal `contigs.csv` rounds to for one
  or two contigs in 20,000. bam-read's per-alignment term reduces exactly to
  `(35 - nm)/35`, so the numerator is carried as two integers and divided
  once (`EXACT_SEQ_TRUE`). Every merged field is now an integer and the
  parallel merge has no rounding to reason about. This is a one-time shift
  against 2.0.0's numbers, at the fifteenth decimal.
- `group_by_fragment` yields decoded `(read, flag, reference_id, length, nm)`
  batches and `score_candidates` takes them, rather than bare
  `AlignedSegment`s (`DECODE_ONCE`). `assign_fragments`, `accumulate_metrics`
  and the CSV writers are unchanged.

### Performance

Profiled, not guessed. Every figure below is from one 1,202,766-record BAM
over 5,000 contigs, and every stage listed produced the identical per-contig
counters — verified by digest, and for the parallel work bit-identical
`p_seq_true` and `p_not_segmented` at 1, 2, 3, 4, 8 and 13 workers.

Per alignment record, through assignment and accumulation — the step that
dominates a real run:

| | µs/record |
| --- | --- |
| 2.0.0 | 2.40 |
| stop paying numpy per alignment and per bin | 1.94 |
| read each record's fields once, not four times | 1.61 |
| decode once across the assign/accumulate boundary | 1.48 |
| slotted accumulator, memoised orphan charge, single-candidate fast path | 1.25 |

1.9× serial, before `--threads` divides what is left across processes.

Per contig, in finalisation — coverage integration, `bases_uncovered` and the
segmenter: **90.3µs → 17.0µs**. Coverage became a difference array integrated
once rather than a slice increment per aligned block; `bin_coverage` stopped
widening the whole coverage vector to int64 per contig and stopped paying
numpy per 30-element bin; the segmenter took a `log` table for its running
state counts and dropped scipy's `logsumexp` for a two-line one; and
`bases_uncovered` counts covered bases rather than building a boolean
temporary to count zeros.

Measured and rejected, recorded so they are not tried again: a `memoryview`
for the coverage bumps (2× faster in isolation, zero end to end), inlining
`_read_length` into `decode` (inside the noise), a `reference_end`-based
CIGAR fast path (~5%, and it changes leading-clip accounting behind a hard
clip), and fusing `iter_alignments` into the grouping loop (cProfile
attributes 0.28µs/record to that generator; it is actually 0.002µs).

`longest_orf` is now the largest per-contig cost at ~57µs — roughly 5s on a
100k-contig assembly, and serial. Two rewrites were tried and neither was
faster. It feeds `n_with_orf`, `mean_orf_percent` and one `contigs.csv`
column, never the score.

### Fixed

- `scripts/compare_orthogroup_picks.py` ordered its `.old.list`/`.new.list`
  by the source it read from — `Orthogroups.txt` line order for
  `--orthogroups`, lexicographic glob order for `--groups` — so the same run
  produced differently ordered lists depending on the flag, and neither
  matched ORP's own `good.<run>.list`. Both now sort by label, which is the
  order ORP keeps deliberately: it reaches contig order in
  `orthomerged.fasta` and so cd-hit-est, where it breaks length ties. Same
  picks, different order.
- `--pick-best` now imports ORP's `best_in_group` as well as `load_scores`
  when the target is ORP 4.0.0 or later, which passes the member list
  directly. Older checkouts keep the mirror.

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

[2.1.0]: https://github.com/macmanes-lab/pytransrate/releases/tag/v2.1.0
[2.0.0]: https://github.com/macmanes-lab/pytransrate/releases/tag/v2.0.0
