# Known issues

Upstream bugs that affect pytransrate runs and that pytransrate cannot work
around. Each entry says how to recognise it, what to do about it now, and what
has to happen before the entry can be deleted.

---

## SNAP-171 — snap 2.0.x SIGFPEs on large runs

| | |
| --- | --- |
| **Upstream** | [amplab/snap#171](https://github.com/amplab/snap/issues/171) |
| **Affects** | every released snap 2.0.x, including the bioconda `snap-aligner=2.0.5` the README pins |
| **Fixed in** | [`0e0997b`](https://github.com/amplab/snap/commit/0e0997b6e2842f7de0be7debdb5acbc275fa7388), released as 2.0.6.dev.2 on snap's `dev` branch — **not in any tag** |
| **Status** | open upstream; one question on the patch unanswered |
| **Workaround** | none available from pytransrate — build snap from `dev` |

### Recognising it

The run dies partway through alignment with:

```
killed by signal 8 (SIGFPE)
```

`logs/snap.log` stops mid-progress-table. The run has usually been going for
some time and the BAM is already large. pytransrate names this issue in the
error text on any SIGFPE; see `SNAP_171` in `src/pytransrate/mapper.py`.

### Cause

A read reaches `computeGlobalScore` with `patternLen == 0`, making
`numVec == 0` at [`AffineGapVectorized.cpp:186`](https://github.com/amplab/snap/blob/3f2017f90207852fda1f73e6ab81cbe82e64dd88/SNAPLib/AffineGapVectorized.cpp#L186)
and dividing by it at [line 351](https://github.com/amplab/snap/blob/3f2017f90207852fda1f73e6ab81cbe82e64dd88/SNAPLib/AffineGapVectorized.cpp#L351).

Bill Bolosky traced the zero in July 2025. When the writer fills its output
buffer partway through writing a read's alignments, it flushes, takes a fresh
buffer and retries the write — and back-clipping from a secondary alignment
was wrongly retained across the retry. A read whose secondary alignment
clipped it to exactly half its length is then clipped to nothing on the retry.
Clipping to any other length produced a silently wrong alignment rather than a
crash.

### Why it looks like a large-assembly bug

It needs three things at once: secondary alignments (`-om`/`-omax`), one of
them back-clipping to exactly half the read, and a buffer refill at that
instant. That makes it a function of **output volume**, not assembly size — a
large assembly is simply what produces a large BAM.

Observed here on a merged assembly of 5,354,958 contigs / 5.66 Gbp, indexed at
`-s 23 -p 1000 -locationSize 5`, 40 threads: snap dies ~17 minutes in with the
BAM already past 100 GB, on every attempt. A full library against that assembly
produces ~160 GB of BAM, where the writer refills continuously and a
per-read-rare coincidence becomes a certainty. The same flags and pipeline on
smaller assemblies run clean — which is why nothing in the test suite catches
it, and why it only appears at the scale where a rerun is expensive.

### Tuning does not avoid it

The original report used the Ruby transrate flags (`-H 300000 -D 5 -om 5`).
pytransrate uses `-H 4000 -D 2 -om 2`, lowered specifically because of an
earlier crash of this kind (see `MULTI_ALIGNMENT_SETTINGS` in `mapper.py`), and
it still dies. Only dropping `-om`/`-omax` entirely helps, and only because it
removes the secondary alignments the bug requires — not an option here, since
those are exactly what `pytransrate.assign` consumes.

Lowering `--multi-edit-distance`, `--max-seed-hits` or
`--max-alignments-per-pair` will not save you.

### What to do

Build snap from the `dev` branch before starting a run that will produce a BAM
in the hundreds of gigabytes:

```bash
git clone -b dev https://github.com/amplab/snap && make -C snap
```

### Open question on the patch

The rewritten loop in `SimpleReadWriter::writePairs` iterates `whichRead` over
both mates but subscripts `clippingForReadAdjustment[0]` for each:

```c
for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
    reads[whichRead]->setAdditionalFrontClipping(result[whichAlignmentPair].clippingForReadAdjustment[0]);
```

2.0.5 used `[0]` for read 0 and `[1]` for read 1, so read 1 now appears to take
read 0's front-clipping adjustment. Asked upstream
([comment](https://github.com/amplab/snap/issues/171#issuecomment-5702151490)),
unanswered as of 2026-09-17. It is a correctness question about clipping, not
about the crash — but until it is settled, treat a `dev` build as the way to
get a large assembly through snap at all rather than as known-good, and sanity
-check its scores against a smaller assembly run on 2.0.5.

### Closing condition

- [ ] Upstream answers the `clippingForReadAdjustment[0]` question
- [ ] snap tags a release containing `0e0997b`
- [ ] bioconda packages that release
- [ ] README install line moves off `snap-aligner=2.0.5`
- [ ] `SNAP_171` note in `mapper.py` and the SIGFPE branch of `_how_it_died`
      updated to name the minimum good version
- [ ] This entry deleted, and the *Known issues* section in `CHANGELOG.md`
      reduced to a line in the release that fixed it
