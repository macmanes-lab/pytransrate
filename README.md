# pytransrate

Quality assessment of de-novo transcriptome assemblies.

pytransrate scores how well an assembly is supported by the reads it was built
from. It maps the reads back, quantifies expression, and reduces the result to
a per-contig score and a single assembly score, so you can compare assemblies
and separate well-supported contigs from junk.

It is a Python port of [transrate](https://github.com/blahah/transrate)
(Smith-Unna et al. 2016), targeting current versions of its dependencies.

## Why the rename

**pytransrate does not reproduce the original's scores, and does not try to.**
The rename exists so nobody compares the two sets of numbers by accident.

The original is pinned to Ruby 2.2.0 (end of life since 2018), ships as an
x86_64-linux-only bundle, and depends on three abandoned binaries. Bringing it
current changed real behaviour:

- **snap-aligner 1.0dev.96 → 2.0.5.** SNAP 2.0 introduced soft clipping, on by
  default. The old `bam-read` helper advanced its reference cursor on soft
  clips as though they consumed reference, which they do not, displacing that
  read's coverage. Inert before, corrupting now. Coverage here follows the SAM
  spec and agrees with `samtools depth`.
- **salmon 0.8.2 → 2.7.0.** `--useErrorModel` is now a hard error, and
  `--sampleOut` is accepted but does nothing, so the `postSample.bam` the old
  pipeline depended on is never written. Fragment assignment moved in-process.
- **Fragment assignment is a redesign, not a port.** salmon 0.8.2 sampled from
  its own posterior; that model is not recoverable from outside. Assignment is
  now a deterministic maximum-a-posteriori choice using salmon's abundances as
  the prior and alignment edit distance as the likelihood.

Several quirks of the published method are preserved deliberately rather than
"fixed", because there is no oracle to say a change is an improvement. They are
documented at named anchors in the source: `BINNING_QUIRK`,
`FRAGMENT_ESTIMATOR`, `CUTOFF_BOUNDARY`, `STATS_QUIRK`, `ORF_CASE_SENSITIVITY`.

## Install

```bash
micromamba create -y -p ./env -c conda-forge -c bioconda \
    python=3.11 numpy scipy pysam snap-aligner=2.0.5 salmon=2.7.0 pip
./env/bin/pip install -e .
```

`conda` or `mamba` work the same way. snap-aligner and salmon must be on PATH.

## Use

```bash
./env/bin/pytransrate -a assembly.fa --left reads.1.fq --right reads.2.fq -t 8 -o results
```

Outputs land in `results/`:

| file | contents |
| --- | --- |
| `assemblies.csv` | one row of assembly-level metrics; `score` and `optimal_score` are columns 37 and 38 |
| `contigs.csv` | per-contig metrics; `score` is column 9 |
| `<assembly>_score_optimisation.csv` | the cutoff/score curve the optimiser walked |

Column *order* in both files is an interface, not presentation — the Oyster
River Protocol reads them positionally. `tests/test_output.py` pins the indices.

### Options worth knowing

| option | why |
| --- | --- |
| `--location-size {4-8}` | snap index sizing. By default pytransrate starts at 4 and steps up on overflow, which costs a full failed index build each time; set it if you already know the value |
| `--seed-size` | the other fix when an assembly overflows the index at every location size |
| `--multi-edit-distance` | snap `-om` (default 2). Controls how many alternate alignments the assignment step gets to choose between |
| `--max-alignments-per-contig` | snap `-mpc` (default 1): best placement per candidate contig |
| `--max-seed-hits` | snap `-H` (default 4000, snap's own default) |
| `--loglevel debug` | logs every external command before it runs |

The snap defaults here differ from the Ruby's, which crash with SIGFPE on real
assemblies; see `MULTI_ALIGNMENT_SETTINGS` in `src/pytransrate/mapper.py`.

`--reference` is **not implemented**. It routed through the unmaintained
`crb-blast` gem. The flag is parsed and raises, rather than silently producing
different output.

## Tests

```bash
./env/bin/python -m pytest
```

The end-to-end tests in `tests/test_pipeline.py` need snap-aligner and salmon
on PATH and skip without them. The rest run anywhere.

Where possible the implementation is checked against something that is not
itself: the segmentation model against a linear-space transcription of the
original C++, coverage against `samtools depth`, and base composition and ORF
length against the original C extension, compiled on demand from
`tests/oracle/orf_oracle.c`.

## Citation

pytransrate implements the method described in the transrate paper. If you use
it, cite the original — see [CITATION.md](CITATION.md).

## License

MIT, as the original. See [LICENSE](LICENSE).
