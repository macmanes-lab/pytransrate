# Changelog

Notable changes to pytransrate. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

Versions up to and including 1.0.3 are the original Ruby
[transrate](https://github.com/blahah/transrate). 2.0.0 is a rewrite in
Python and **does not reproduce their scores** — see *Changed* below.

## [Unreleased]

### Fixed

- **The read-metrics step no longer asks for more memory than the machine
  has.** It divides the fragments across `--threads` processes, and each one
  holds its own full-size copy of the per-base coverage accumulators — 4
  bytes a base. On a 5.8 Gbp merged assembly that is 23 GB per process, so
  `-t 40` asked for 928 GB and the OOM killer ended the run *after* nine
  hours of mapping and quantifying:

  ```
  assigning across 40 processes (927.7 GB of shared accumulators)
  Killed
  ```

  The budget is now worked out first — from the cgroup, the Slurm
  allocation and `/proc/meminfo`, or from the new `--max-memory`/`--mem` —
  and the processes are capped at what fits, with a warning naming both the
  figure and where it came from:

  ```
  40 processes would need 927.7 GB of shared accumulators and the Slurm
  allocation is 773.1 GB; assigning across 25 instead.
  ``` snap and
  salmon still use every thread. Little speed is lost: dividing this step is
  bounded at about 6x however many processes run (see STRIDING in
  `read_metrics`), which is reached around 16.

  The parent also frees each worker's buffer as soon as it has been summed
  rather than at the end of the step, so its own copy is taken while two
  buffers are held instead of `--threads` plus one.

- **A killed run no longer half-reuses what it left behind.** Resuming
  already reused the snap index, the BAM and `salmon/quant.sf`, which is
  what makes a retry cheap — but a process killed mid-write leaves a BAM
  with no BGZF end-of-file marker and a `quant.sf` missing contigs, and both
  were reused without a word, giving metrics computed against a fraction of
  the library. Both are now checked before they are trusted: the BAM for its
  end-of-file marker, `quant.sf` for a complete final row and one row per
  contig in the BAM. Anything short is redone.

  Every reuse decision is also logged now, the way the index build has been
  since it was worth knowing which of them a run skipped. Reusing the BAM
  silently was the difference between a five-hour step and no step at all,
  with nothing in the log either way.

### Added

- **Notes on all of this**, since a memory ceiling nobody can see coming is
  the kind of thing that gets rediscovered: `README.md` gains the per-process
  figure and the resume behaviour, `USAGE.md` covers both ways to set the
  budget and what the warning means, MEMORY_BUDGET in `read_metrics` records
  why the accumulators are per-worker at all, and a test pins the warning
  USAGE quotes to the one the code prints.

- **`--max-memory SIZE`, also spelled `--mem SIZE`.** What the read-metrics
  step may use — `670G`, `512M`, or a bare number for GB (`670Gi` for GiB).
  Only needed where the detected figure is wrong, which is most likely on a
  scheduler enforcing a limit the cgroup does not show, or under a pipeline
  that already knows what it asked for: `--mem` is what those call it, so
  the figure can be passed straight through. An explicit figure wins over
  detection, and reads back unchanged in the log, decimal as everything else
  this program prints is.

- **Gzipped input.** The assembly and the read files may each be gzipped, in
  any combination: `pytransrate -a asm.fa.gz --left r1.fq.gz --right
  r2.fq.gz`. Compression is detected by the magic number rather than the file
  name, so a gzipped file named anything is read correctly and a plain file
  named `.gz` is not mangled; bgzipped input works for the same reason.

  The read files are handed to snap-aligner still compressed — it reads
  gzipped FASTQ itself, and a library is routinely tens of gigabytes. The
  assembly cannot be: snap and salmon take a filename, not a handle, so a
  gzipped assembly is decompressed into the output directory for the length
  of the index-and-quantify step and removed afterwards. A run without reads
  decompresses nothing at all, since the FASTA is read in-process.

  Scores are unaffected — the pipeline test asserts a gzipped run reproduces
  an uncompressed one contig for contig.

- **`logs/snap.log` is now written as snap runs**, and covers indexing as
  well as mapping. It was assembled in memory and written once snap returned,
  so the runs whose log is actually wanted — snap killed by an OOM killer or
  a scheduler wall clock, or a `^C` — left no log at all. snap now writes
  straight to the file, which survives anything short of the machine going
  away. Each entry is preceded by the command that produced it, so a sweep
  across several `-locationSize` values stays readable, and errors raised out
  of the mapping step quote the tail of the output and name the log rather
  than reproducing snap's whole progress table.

  A snap that dies on a signal now says so by name — `killed by signal 8
  (SIGFPE)` rather than `exit -8` — and says which way to look: a crash
  points at the multiple-alignment flags (see MULTI_ALIGNMENT_SETTINGS),
  whereas a SIGKILL or SIGTERM points at the job's memory ceiling or wall
  clock instead.

### Fixed

- **A crashing snap no longer takes its own last words with it.** Writing the
  log as snap ran was only half the problem: the buffering that swallowed the
  output is in snap, not in pytransrate. C stdio block-buffers a few kilobytes
  when its output is a file rather than a terminal and flushes at exit, and a
  process killed by a signal never exits — so a snap-aligner 2.0.5 run that
  took SIGFPE twenty minutes into a merged assembly left a log containing one
  line, `Welcome to SNAP version 2.0.5.`, which survived only because snap
  writes its banner to stderr. The index load, the bases indexed and the
  progress table were all still sitting in the buffer, which is to say there
  was no way to tell whether snap died loading the index or partway through
  the reads.

  snap is now run under `stdbuf -oL -eL` where coreutils provides it — every
  Linux cluster, not macOS — so each line reaches the log as it is printed and
  a crash leaves the position it crashed at. Where `stdbuf` is missing the run
  proceeds unchanged, with the log no more truncated than it already was.

- **Contig names are no longer cut at the first `|`.** They are taken from
  the defline up to the first space, which is the rule snap-aligner and
  salmon apply when they take a reference name from a FASTA.

  The Ruby cut at `|` as well, following BioRuby's `entry_id`, and so did
  this port. That made any assembly with pipes in its deflines unusable
  twice over. ENA and TSA downloads — `>ENA|GADU01000001|GADU01000001.1
  Gadus morhua mRNA` — collapse to the single identifier `ENA` and were
  rejected with `Non-unique fasta identifier found: >ENA` before any work
  started. Worse, an assembly whose truncated names did happen to stay
  unique ran to completion and matched nothing in the BAM header or in
  salmon's `Name` column, so every read metric in `contigs.csv` came back
  zero with no error anywhere.

  Assemblies without pipes — Trinity, rnaSPAdes, the ORP's merged output —
  are unaffected: their names never contained one to cut at. The duplicate
  check and its message remain, for deflines that genuinely collide in
  their first field. The pipeline test now runs the same dataset a second
  time under ENA-style deflines and asserts the scores match the plain-named
  run contig for contig, so the join is checked against the binaries rather
  than taken on trust.

- **A `-locationSize` failure now reports what is actually filling the
  genome.** snap sizes a genome as `fileSize + (nContigs + 1) * padding` and
  applies the ceiling to that total, so on a fragmented transcriptome most of
  what overflows is padding rather than assembly — at 5.4M contigs the default
  padding alone is 10.7 Gbp against a 4-byte ceiling of 4.29 Gbp. The warning
  now gives the FASTA size, the contig count, the padding's share of the
  total, and the ceiling, so `--padding` is visible as the lever it is.
  Computed once per build and only on the failure path.

- **A completed index is never deleted by the `-locationSize` sweep.** The
  sweep clears the directory before retrying at a larger size, which is right
  for the partial build it just made and wrong for a finished index that was
  already there. It now refuses and says so instead. Reuse and rebuild are
  also logged at INFO with the resolved path, because they were
  indistinguishable in a log: the reuse branch returned silently, so a run
  that rebuilt an index it should have reused looked exactly like a first run.

- **An index build that exits 0 without writing an index is now caught.** snap
  can do this — it is the same silent failure the mapping step already guards
  against — and mapping against an index that is not there dies with nothing
  but snap's version banner, which is very hard to read backwards. The build
  checks for `GenomeIndex` before reporting success.

- **An index is no longer deleted while another run is aligning against it.**
  The index directory is named after the assembly, so two runs of the same
  assembly into the same output directory shared it; if the second reached
  the `-locationSize` retry — the one place an index is ever deleted — while
  the first was already mapping, snap lost its genome mid-alignment. An
  exclusive `flock` is now held on `<assembly>.index.lock` beside the
  directory, from the start of the index build until the run is done mapping,
  and the second run stops with a message naming the lock instead. The kernel
  drops the lock when the holder dies, so a killed run leaves nothing to
  clean up by hand.

- The `-locationSize` sweep now runs on every way snap reports that the
  location size is too small. It matched two of snap's four messages, so a
  large assembly failing with `Genome is too big for 4 byte genome
  locations.  Specify a larger location size with -locationSize` — the check
  snap makes before any index work, and the one a big merged assembly hits
  first — was reported as a hard failure at `-locationSize 4` instead of
  being retried at 5. `Not enough address space to index this genome with
  this seed size` was missed the same way; despite its wording it is the
  same location-namespace limit, bounded by `2**(locationSize*8) - 1`, not
  a machine memory limit.

- The read-count fallback counted lines of compressed data when a BAM was
  reused and the saved count file was missing, giving a meaningless number of
  fragments for gzipped reads. It decompresses first. The path only runs when
  `--keep-bam` output is reused without its `*-read_count.txt`.

### Changed

- **snap's contig padding drops from 2000 to 1000**, with `--padding` added to
  override it. snap pads every contig with Ns so an alignment cannot run off
  one contig into the next, and counts those Ns as genome: `FASTA.cpp` sizes
  the genome as `fileSize + (nContigs + 1) * padding`, and the `-locationSize`
  ceiling of `2**32 - 16` bases applies to that total.

  snap's default of 2000 is sized for genomes, where a few hundred contigs
  make the padding a rounding error. A transcriptome inverts it: at 1.5M
  contigs averaging a kilobase, the padding contributes 3 Gbp of Ns on top of
  ~1.5 Gbp of sequence, so the padding is larger than the assembly and is on
  its own enough to cross the 4-byte ceiling. Crossing it forces the
  `-locationSize` sweep up to 5, which pays for a second full index build and
  then leaves the aligner holding a larger index in memory for the run.

  1000 is the largest reduction that costs nothing on either bound snap
  documents for this value. The correctness floor is the maximum edit
  distance, which is 30 — two orders of magnitude clear. The other is a
  performance note about the padding exceeding the paired-end gap, which is
  the `-s` maximum of 1000: padding of 1000 sits at that bound rather than
  above it, so a pair straddling two adjacent contigs is no longer separated
  by more than the maximum spacing. It cannot be called a proper pair either
  way — crossing the padding means crossing 1000 Ns, which no alignment
  within an edit distance of 30 survives — so what changes is snap doing the
  work to reject it, not the rejection. Pass `--padding` above 1000 + read
  length if that ever shows up in a profile.

  **This moves scores on any assembly where it changes the index**, since
  contigs land at different genome locations. An index already built at 2000
  is still read as-is; the padding is stored in it, and `build_index` reuses
  a complete index rather than rebuilding it.

- The banner is blue and yellow, replacing the green/yellow/red flanks
  inherited from the Ruby. Colour remains opt-in and off when `NO_COLOR` is
  set, `TERM=dumb`, or stderr is not a terminal.

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
