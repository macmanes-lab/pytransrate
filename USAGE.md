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

### Gzipped input

The assembly and the read files may be gzipped, in any combination, and may
be mixed with uncompressed ones:

```bash
pytransrate -a assembly.fa.gz --left r1.fq.gz --right r2.fq.gz -o results
```

Compression is detected by reading the first two bytes, not by the file
name, so a gzipped file called something other than `.gz` is read correctly
and a plain file called `.gz` is not mangled. bgzipped input works too, being
valid gzip.

The read files are passed to snap-aligner compressed, which is why a library
of any size costs nothing here. The assembly is different: snap and salmon
are given a filename rather than a handle, so a gzipped assembly is
decompressed once into the output directory for them and deleted when they
are done — briefly, one extra copy of the assembly on disk. A run without
reads decompresses nothing, since pytransrate reads the FASTA itself.

Output naming ignores the suffix: `asm.fa.gz` gives the same result
directory and the same `asm.fa_score_optimisation.csv` as `asm.fa` would.
The `assembly` column of `assemblies.csv` records the path as given.

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
| `-t`, `--threads N` | 8 | threads for snap and salmon, and processes for the read-metrics step |
| `--max-memory SIZE`, `--mem SIZE` | detected | memory the read-metrics step may use — `670G`, `512M`, or a bare number for GB (`670Gi` for GiB). It caps the processes that step forks, since each holds a full-size copy of the coverage accumulators. The default reads the cgroup, the Slurm allocation and `/proc/meminfo`, so the flag is only needed when that figure is wrong |
| `--loglevel LEVEL` | `info` | `error`, `warn`, `info`, `debug`; `debug` logs every external command before it runs |
| `--keep-bam` | off | keep the alignment BAM instead of deleting it on success. It is large — many gigabytes on a real library |
| `--no-banner` | off | suppress the startup banner |
| `--version` | | print the version and exit |

### snap index

Only needed when a large or repetitive assembly exhausts snap's genome
locations.

| option | default | |
| --- | --- | --- |
| `--location-size {4-8}` | sweep 4→8 | bytes per genome location. The default retries upward whenever snap runs out of locations, and each failed attempt is a full index build — set this if you already know the value |
| `--seed-size N` | 23 | the other fix when an assembly is still too big at every location size |
| `--padding N` | 1000 | snap `-p`, Ns inserted between contigs. Counted toward the genome size, so on an assembly with millions of contigs the padding, not the sequence, is usually what exhausts the location namespace. Must stay above `--edit-distance` |

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

**The index runs out of genome locations.** snap says one of four things:

- `Genome is too big for 4 byte genome locations` — the assembly itself is
  larger than the location namespace. Big merged assemblies hit this one.
- `Ran out of overflow table namespace`
- `Trying to use too many overflow entries`
- `Not enough address space to index this genome with this seed size` —
  despite the wording this is the same namespace limit, not a RAM limit.

All four mean the same thing, and pytransrate retries at `-locationSize`
4, 5, 6, 7, 8 in turn on any of them. Each failed attempt is a complete index
build, so if you already know an assembly needs 6, pass `--location-size 6`
and skip the wasted work. If every size is still too small, raise
`--seed-size`.

**Check how much of that "genome" is padding first.** snap pads every contig
with Ns and counts them as genome, so the size the location namespace has to
cover is `n_bases + padding * (n_contigs + 1)`, not the assembly. On a
transcriptome with millions of short contigs the padding is routinely larger
than the sequence:

```
grep -c '^>' assembly.fa                          # n_contigs
grep -v '^>' assembly.fa | tr -d '\n' | wc -c     # n_bases
```

If `n_bases` alone is comfortably under 4,294,967,280 but `n_bases + 1000 *
n_contigs` is not, lowering `--padding` is the cheaper fix than stepping
`--location-size` up: it avoids a second full index build and leaves a
smaller index for the aligner to hold in memory. The floor is
`--edit-distance` (30 by default). Note that changing the padding changes
where contigs sit in the genome, and so moves scores.

**A very redundant assembly.** Duplicate contigs differing by 0–2 bases
generate large numbers of secondary alignments. The defaults are already
tuned for this — `-om 2` rather than the Ruby's 5 retains 98.6% of secondary
alignments and 99.98% of fragments hitting more than one contig, measured on
a redundant 27,000-contig assembly. Raising `--multi-edit-distance` costs a
great deal of time for very little extra signal.

**Threads, and the memory they cost.** snap and salmon both take `-t`, and
for them memory scales with the index, not the thread count. The read-metrics
step is different: it divides the fragments across processes, and each one
needs its own copy of the per-base coverage accumulators, at 4 bytes a base:

```
bytes per process ~= 4 * (n_bases + n_contigs)
```

A 5.8 Gbp merged assembly is 23 GB per process, so `-t 40` asks for 928 GB —
and since mapping and quantifying come first, the OOM killer arrives hours
into the run. pytransrate now works out what is available (cgroup, Slurm
allocation, `/proc/meminfo`), caps the processes at what fits, and says so,
naming the figure it used and where that figure came from:

```
[ WARN] 40 processes would need 927.7 GB of shared accumulators and the Slurm
        allocation is 773.1 GB; assigning across 25 instead. Pass --max-memory
        (or --mem) to say otherwise; mapping and quantifying still used every
        thread.
```

snap and salmon still get every thread; only this step is capped. Little is
lost — dividing this step is bounded at about 6x however many processes run,
so it is at its plateau by ~16.

Two ways to set the budget, and both matter:

- **On its own**, pytransrate detects it: cgroup limit, then the Slurm
  allocation (`SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` × CPUs), then
  `MemAvailable`. It takes the smallest, since every one of them is real, and
  names the one it used — `#SBATCH --mem 720G` reads back as `the Slurm
  allocation is 773.1 GB`, which is the same figure counted in GB rather than
  GiB. Where that name is `available memory (MemAvailable)` on a shared node,
  the figure is whatever the node had free at that moment and is worth
  overriding.
- **From a pipeline that already knows the figure**, pass it through:
  `--max-memory 670G`, or `--mem 670` for the same thing, since that is what
  the wrapper and the scheduler directive call it. An explicit figure wins
  over detection — which is the point, because the case that needs it is a
  scheduler enforcing a limit the cgroup does not show.

Sizes are decimal unless you ask otherwise: `670G` is 670 GB, `670Gi` is 670
GiB, and a bare `670` is 670 GB. The figure you pass is the figure the log
reports back.

**Resuming a killed run.** A run that dies after mapping does not repeat it.
Rerunning the same command with the same `-o` reuses, in order:

- the snap index, if its directory holds `GenomeIndex`;
- the BAM, if it ends with the BGZF end-of-file marker — a BAM left by a
  killed snap does not, and is mapped again rather than half-read;
- `salmon/quant.sf`, if it ends on a complete row and carries one row per
  contig in the BAM.

Each decision is logged either way, so the log says which hours were skipped
and which were spent. The BAM is deleted on a *successful* run unless
`--keep-bam` is given; on a failed one it is left exactly so the next attempt
can use it.

**Very large assemblies need a patched snap.** Anything producing a BAM in the
hundreds of gigabytes will hit [amplab/snap#171][snap171] and die with SIGFPE
partway through alignment, whatever the flags. Build snap from its `dev` branch
(2.0.6.dev.2 or later) before starting a run of that size; see
*Troubleshooting* below.

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

**Killed, exit 137, right after `assigning fragments and computing read
metrics`** — the OOM killer. The line above it reports what that step asked
for: `assigning across 40 processes (927.7 GB of shared accumulators)`.
Current versions cap the processes at what memory allows instead, and the
warning they print names both the figure and where it came from. If that
figure is wrong for your scheduler, pass `--max-memory 670G` (or `--mem 670`,
the same option). Rerun the same command with the same `-o` — the index, BAM
and `quant.sf` are all reused, so a retry costs minutes rather than the hours
the first attempt spent before it died.

**A run failed — do I have to map again?** No. Rerun the same command with
the same `-o`: the index, the BAM and `quant.sf` are all reused. The BAM is
only deleted once a run has succeeded all the way to `assemblies.csv`, so a
failed run keeps it. A reusable BAM has an `<bam>.align.done` beside it naming
the command that wrote it.

**`<bam>.partial` in my output directory** — a previous run was killed partway
through mapping. It is kept rather than overwritten, because it is the evidence
of what the aligner did before it died, but it cannot be used for metrics and
nothing will read it. Delete it whenever you like; at most one is kept.

**snap dies with SIGFPE** — this is [amplab/snap#171][snap171], a
divide-by-zero in snap 2.0.x, and **the defaults do not avoid it**. Upstream
fixed it in 2.0.6.dev.2; until that reaches a release you need snap built from
the `dev` branch:

```bash
git clone -b dev https://github.com/amplab/snap && make -C snap
```

The bug needs secondary alignments (`-om`/`-omax`) and the writer running out
of output buffer in the middle of writing a read, so it tracks how big the BAM
is, not how big the assembly is — it just takes a large assembly to produce a
large BAM. Past roughly 100 GB of output it stops being unlucky and becomes
reliable. Lowering `--multi-edit-distance`, `--max-seed-hits` or
`--max-alignments-per-pair` will not save you; only dropping `-om`/`-omax`
altogether does, and pytransrate cannot do that because fragment assignment is
what consumes those secondary alignments. See [the changelog][known] for the
full account, including one open question about the upstream patch.

[snap171]: https://github.com/amplab/snap/issues/171
[known]: CHANGELOG.md

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
