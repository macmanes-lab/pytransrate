# pytransrate usage

Full reference. For what the tool is and how to install it, see
[README.md](README.md); for what changed from the Ruby transrate, see
[CHANGELOG.md](CHANGELOG.md).

- [Running](#running)
- [Output](#output)
- [What the metrics mean](#what-the-metrics-mean)
- [Options](#options)
- [Tuning for hard assemblies](#tuning-for-hard-assemblies)
- [Comparing runs](#comparing-runs)
- [Troubleshooting](#troubleshooting)
- [Use from the Oyster River Protocol](#use-from-the-oyster-river-protocol)

## Running

```bash
pytransrate -a assembly.fa --left r1.fq --right r2.fq -t 16 -o results
```

`--left` and `--right` must be given together and matched in order. Multiple
files are comma-separated lists of the same length:

```bash
pytransrate -a assembly.fa --left a1.fq,b1.fq --right a2.fq,b2.fq -o results
```

Several assemblies in one invocation each get their own subdirectory of the
output, and share one `assemblies.csv` with a row apiece:

```bash
pytransrate -a one.fa,two.fa --left r1.fq --right r2.fq -o results
```

Without reads only sequence metrics are computed — no aligner or quantifier
is needed, and the read columns are absent from both CSVs:

```bash
pytransrate -a assembly.fa -o results
```

The output directory is refused if it already holds an `assemblies.csv`, so a
rerun cannot silently overwrite a previous result. Give each run its own `-o`.

### The run report

Everything goes to **stdout**; the banner goes to stderr. The first line
records the invocation, shell-quoted, so a log identifies the settings that
produced it:

```
[ INFO] 2026-09-07 09:33:47 : command: pytransrate -a assembly.fa --left r1.fq --right r2.fq -o results
[ INFO] 2026-09-07 09:33:47 : loading assembly: assembly.fa
[ INFO] 2026-09-07 09:33:47 : calculating contig metrics...
[ INFO] 2026-09-07 09:33:47 : contig metrics:
[ INFO] 2026-09-07 09:33:47 :   min contig length  201
[ INFO] 2026-09-07 09:33:47 :   max contig length  12,394
[ INFO] 2026-09-07 09:33:47 :   N50                1,874
```

then, once reads are mapped, the complete read-metrics block, the soft-clip
summary, and the scores. Numbers are formatted exactly as `assemblies.csv`
rounds them, so a value read off the log is the value in the file.

## Output

### `assemblies.csv` — one row per assembly, 40 columns

| # | column | |
| --- | --- | --- |
| 1 | `assembly` | path as given |
| 2–19 | `n_seqs` `smallest` `largest` `n_bases` `mean_len` `n_under_200` `n_over_1k` `n_over_10k` `n_with_orf` `mean_orf_percent` `n90` `n70` `n50` `n30` `n10` `gc` `bases_n` `proportion_n` | sequence only; computed from the FASTA, no reads involved |
| 20–36 | `fragments` `fragments_mapped` `p_fragments_mapped` `good_mappings` `p_good_mapping` `bad_mappings` `potential_bridges` `bases_uncovered` `p_bases_uncovered` `contigs_uncovbase` `p_contigs_uncovbase` `contigs_uncovered` `p_contigs_uncovered` `contigs_lowcovered` `p_contigs_lowcovered` `contigs_segmented` `p_contigs_segmented` | read mapping |
| 37 | `score` | the assembly score |
| 38 | `optimal_score` | best score reachable by discarding low-scoring contigs |
| 39 | `cutoff` | contig score at which to cut to reach `optimal_score` |
| 40 | `weighted` | expression-weighted variant of `score` |

Floats are rounded to 5 places.

### `contigs.csv` — one row per contig, 18 columns

| # | column | |
| --- | --- | --- |
| 1–4 | `contig_name` `length` `prop_gc` `orf_length` | sequence only |
| 5–14 | `in_bridges` `p_good` `p_bases_covered` `p_seq_true` `score` `p_not_segmented` `eff_length` `eff_count` `tpm` `coverage` | read mapping |
| 15–18 | `sCnuc` `sCcov` `sCord` `sCseg` | the four score components |

Floats are rounded to 6 places. Columns 15–18 duplicate 8, 7, 6 and 10 under
the names used in the paper; both sets are written, as the Ruby did.

**The column order is an interface.** Downstream tools read these files
positionally — ORP takes `score`/`optimal_score` from zero-based indices
36/37 of `assemblies.csv` (columns 37 and 38 above) and the contig score from
column 9 of `contigs.csv`.
`tests/test_output.py` pins them. Do not reorder.

### `<assembly>_score_optimisation.csv`

The cutoff/score curve the optimiser walked: one row per candidate cutoff,
with the assembly score that cutting there would produce.

## What the metrics mean

### The contig score

```
score = sCnuc × sCcov × sCord × sCseg      each floored at 0.01
```

| component | column | what it measures | it falls when |
| --- | --- | --- | --- |
| `sCnuc` | `p_seq_true` | per-base accuracy, from alignment edit distance | reads disagree with the contig's sequence |
| `sCcov` | `p_bases_covered` | fraction of bases with at least one read over them | parts of the contig have no read support |
| `sCord` | `p_good` | fraction of fragments mapping consistently | mates land wrongly oriented, too far apart, or on a different contig |
| `sCseg` | `p_not_segmented` | probability the coverage is one transcript, not several | coverage steps up or down mid-contig, as when two transcripts are joined |

Because it is a product, one bad component sinks a contig regardless of the
others. That is the intent: a contig needs to be right in every respect.

### The assembly score

```
score = geomean(contig scores) × (good_mappings / fragments)
```

Two independent things: how good the contigs are, and how much of the library
the assembly explains. A score can move because either changed, and they mean
different things — `scripts/compare_transrate_runs.py` splits any difference
into these two exactly.

`optimal_score` walks contigs worst-first, removing each in turn, and records
the assembly score that would result. The best is `optimal_score`, and
`cutoff` is the contig score at which to cut. A large gap between `score` and
`optimal_score` means the assembly carries a lot of junk that filtering would
remove.

Scores are comparable **between assemblies built from the same reads**. They
are not comparable across libraries — a score depends on read depth and read
quality as much as on assembly quality — and not comparable to scores from the
Ruby transrate.

### Assembly-level read metrics

| column | |
| --- | --- |
| `fragments` | fragments in the library |
| `fragments_mapped` / `p_fragments_mapped` | fragments with at least one alignment |
| `good_mappings` / `p_good_mapping` | fragments mapping consistently: both mates on the same contig, correctly oriented, plausibly spaced |
| `bad_mappings` | mapped but not good |
| `potential_bridges` | fragments whose mates landed on different contigs |
| `bases_uncovered` / `p_bases_uncovered` | assembly bases with no read over them |
| `contigs_uncovbase` | contigs with at least one uncovered base |
| `contigs_uncovered` | contigs whose mean coverage is under 1× |
| `contigs_lowcovered` | contigs whose mean coverage is under 10× — this *includes* the uncovered ones |
| `contigs_segmented` | contigs whose coverage looks like more than one transcript |

## Options

### Required

| option | |
| --- | --- |
| `-a`, `--assembly FASTA` | assembly file(s), comma-separated |

### Reads

| option | |
| --- | --- |
| `--left FASTQ` | left reads, comma-separated |
| `--right FASTQ` | right reads, comma-separated |

Give both or neither.

### General

| option | default | |
| --- | --- | --- |
| `-o`, `--output DIR` | `transrate_results` | output directory |
| `-t`, `--threads N` | 8 | threads for snap and salmon |
| `--loglevel LEVEL` | `info` | `error`, `warn`, `info`, `debug`; `debug` logs every external command before it runs |
| `--keep-bam` | off | keep the alignment BAM instead of deleting it on success. It is large — many gigabytes on a real library |
| `--no-banner` | off | suppress the startup banner |
| `--version` | | print the version and exit |

### snap index

Only needed when a large or repetitive assembly overflows the index.

| option | default | |
| --- | --- | --- |
| `--location-size {4-8}` | sweep 4→8 | bytes per genome location. The default retries upward on overflow, and each failed attempt is a full index build — set this if you already know the value |
| `--seed-size N` | 23 | the other fix when an assembly overflows at every location size |

### snap mapping

Defaults are the configuration verified against real assemblies. The Ruby's
values crash snap 2.x with SIGFPE.

| option | snap flag | default | |
| --- | --- | --- | --- |
| `--max-alignments-per-contig` | `-mpc` | 1 | best placement per candidate contig, applied before `-omax`. **Do not raise it** — see below |
| `--max-alignments-per-pair` | `-omax` | 10 | cap on alignments per pair |
| `--multi-edit-distance` | `-om` | 2 | extra edit distance admitted for secondary alignments, which is what fragment assignment chooses between |
| `--extra-search-depth` | `-D` | 2 | must be ≥ `--multi-edit-distance` |
| `--max-seed-hits` | `-H` | 4000 | snap's own default; the Ruby used 300000 with no recorded rationale |
| `--edit-distance` | `-d` | 30 | max edit distance per pair |
| `--max-candidate-pool` | `-mcp` | not passed | must be under 2147483647; the Ruby's value overflowed snap's `atoi()` |

**`--max-alignments-per-contig` raises the score without improving the
assembly.** Measured on three assemblies of one library, disabling or
loosening the cap raises the transrate score by 0.003–0.005, and 98–105% of
that comes from `fragments_mapped` and `good_mappings` counting the same
fragment more than once on the contig it was assigned to. `sCnuc` falls at the
same time. The accumulator documents one alignment per fragment per reference
as its precondition, and `-mpc 1` is what guarantees it. Full measurement at
`MULTI_ALIGNMENT_SETTINGS` in `src/pytransrate/mapper.py`.

### salmon

| option | |
| --- | --- |
| `--no-error-model` | don't pass `--errorModel`. Only sensible if your BAM carries `AS` tags; snap-aligner does not emit them |

### Not implemented

`-r`, `--reference FASTA` — reference-based (CRB-BLAST) metrics. The Ruby
routed these through the unmaintained `crb-blast` gem. The flag is parsed and
raises, so a script using it gets a clear error rather than silently different
output.

## Tuning for hard assemblies

**The index overflows.** snap reports `Ran out of overflow table namespace` or
`Trying to use too many overflow entries`. pytransrate retries at
`-locationSize` 4, 5, 6, 7, 8 in turn, and each failed attempt is a complete
index build. If you already know an assembly needs 6, pass
`--location-size 6` and skip the wasted work. If every size overflows, raise
`--seed-size`.

**A very redundant assembly.** Duplicate contigs differing by 0–2 bases
generate large numbers of secondary alignments. The defaults are already
tuned for this — `-om 2` rather than the Ruby's 5 retains 98.6% of secondary
alignments and 99.98% of fragments hitting more than one contig, measured on
a redundant 27,000-contig assembly. Raising `--multi-edit-distance` costs a
great deal of time for very little extra signal.

**Threads.** snap and salmon both take `-t`. Memory scales with the index, not
the thread count.

## Comparing runs

Two scripts, both reading only the CSVs, so they cannot perturb what they
measure.

### `scripts/compare_transrate_runs.py`

```bash
python scripts/compare_transrate_runs.py base=run_a alt=run_b --baseline base
```

Prints five sections:

1. **Validity** — the sequence-only columns must match, since they come from
   the FASTA alone. If they don't, the runs analysed different assemblies and
   nothing else is meaningful. Compared numerically, so two implementations
   disagreeing in the last decimal is reported as such rather than as a
   failure.
2. **Assembly-level metrics**, with deltas against the baseline.
3. **Score decomposition** — splits each score difference exactly into its
   good-rate and contig-geomean causes.
4. **Per-contig paired deltas** for `score` and its four components, matched
   on contig name.
5. **Rank agreement** — Pearson *r*, the share of random contig pairs the two
   runs order oppositely, and whether each run's own cutoff keeps the same
   contigs.

Pass `--replicate LABEL` naming a rerun at the baseline's settings and every
difference is expressed as a multiple of the run-to-run noise. Without one
there is no noise floor; snap and salmon are both multithreaded and neither
promises identical output run to run.

Each comparison also emits one `SUMMARY` line, so a sweep over datasets
reduces to `| grep ^SUMMARY`.

### `scripts/compare_orthogroup_picks.py`

```bash
python scripts/compare_orthogroup_picks.py old/contigs.csv new/contigs.csv \
    --orthogroups /path/to/Orthogroups.txt --out-prefix picks
```

Runs ORP's orthogroup selection against both runs and reports exactly which
groups change representative, split by group size and by how decisively each
run preferred its own winner. Use the `contigs.csv` from the run over
`merged.fasta`, not over a finished assembly.

## Troubleshooting

**`the aligner wrote a BAM record htslib cannot parse`** — you are on an old
build. Current versions skip such records, count them, and report
`N alignment(s) skipped` at WARNING. It is a snap bug at contig boundaries,
present in every 2.0.x release, and the count is normally in the single
digits.

**`100 unreadable records in a row`** — that is a truncated or corrupt BAM,
usually an aligner run that died partway. Delete it and map again.

**snap dies with SIGFPE** — you have passed the Ruby's multi-alignment
settings. Use the defaults.

**`assemblies.csv would be overwritten`** — give the run its own `-o`.

**Which build am I running?**

```bash
pytransrate --version
```

**Scores don't match a previous run** — check the first line of both logs.
Every run records its own invocation.

## Use from the Oyster River Protocol

ORP invokes transrate with no mapping flags, so every default above applies,
and reads the results positionally: `score` and `optimal_score` from zero-based
indices 36 and 37 of `assemblies.csv`, and the contig score from column 9 of
`contigs.csv` in `scripts/pick_best_contigs.py`, which keeps the
highest-scoring member of each orthogroup.

That last point matters when changing anything that affects contig scores: it
is the *ordering* of contig scores that reaches the assembly, not their level.
Two configurations can differ substantially in score and select identically,
or agree closely and not. `scripts/compare_orthogroup_picks.py` measures it
directly.
