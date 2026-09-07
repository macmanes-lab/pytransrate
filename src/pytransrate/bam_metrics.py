"""Per-contig metrics from a BAM alignment.

Replaces the external ``bam-read`` binary (Blahah/transrate-tools), which the
Ruby transrate shelled out to.  Produces the same ten columns:

    name, p_seq_true, bridges, length, fragments_mapped, both_mapped,
    properpair, good, bases_uncovered, p_not_segmented

One deliberate behavioural change from the C++ is documented at SOFT_CLIP_FIX.
Two quirks are preserved on purpose and documented at FRAGMENT_ESTIMATOR and
in :mod:`~pytransrate.segmenter` (BINNING_QUIRK).

The input is expected to carry exactly one alignment per fragment -- the role
``postSample.bam`` played under salmon 0.8.2's ``--sampleOut``.  See
:mod:`~pytransrate.assign`.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np
import pysam

from pytransrate.segmenter import DEFAULT_NULL_PRIOR, bin_coverage, prob_not_segmented

logger = logging.getLogger("pytransrate")

__all__ = [
    "CSV_COLUMNS",
    "MalformedBamError",
    "MalformedRecordStats",
    "iter_alignments",
    "ContigMetrics",
    "estimate_realistic_distance",
    "accumulate_metrics",
    "compute_bam_metrics",
    "write_metrics_csv",
]

#: Output columns, in the order bam-read emitted them.
CSV_COLUMNS = (
    "name",
    "p_seq_true",
    "bridges",
    "length",
    "fragments_mapped",
    "both_mapped",
    "properpair",
    "good",
    "bases_uncovered",
    "p_not_segmented",
)

# CIGAR opcodes as pysam reports them in ``cigartuples``.
_CIGAR_MATCH = 0  # M
_CIGAR_INS = 1  # I
_CIGAR_DEL = 2  # D
_CIGAR_REF_SKIP = 3  # N
_CIGAR_SOFT_CLIP = 4  # S
_CIGAR_HARD_CLIP = 5  # H
_CIGAR_PAD = 6  # P
_CIGAR_EQUAL = 7  # =
_CIGAR_DIFF = 8  # X

#: Ops that consume reference *and* represent aligned bases.
_COVERING_OPS = frozenset({_CIGAR_MATCH, _CIGAR_EQUAL, _CIGAR_DIFF})

#: Ops that advance the reference cursor without contributing coverage.
_SKIPPING_OPS = frozenset({_CIGAR_DEL, _CIGAR_REF_SKIP})

# ---------------------------------------------------------------------------
# SOFT_CLIP_FIX
#
# pileup.cpp advanced the reference cursor on soft clips:
#
#     if (op.Type == 'D' || op.Type == 'S') {
#       pos += (int)op.Length;
#     }
#
# Per the SAM spec, S consumes *query*, not reference -- and BAM's POS is
# already the leftmost *mapped* base, so the clip is excluded before we start.
# Advancing on S therefore shifts a read's whole coverage contribution
# rightward by the clip length.
#
# This was inert under snap-aligner 1.0dev.96, which never emitted soft clips.
# SNAP 2.0.0 introduced them, so preserving the C++ behaviour would silently
# corrupt bases_uncovered (sCcov) and the binned states (sCseg) -- both direct
# multiplicands of the contig score.
#
# We follow the spec.  Unlike the segmenter's binning quirk this is a
# demonstrable defect with an external oracle: `samtools depth -a` agrees with
# the handling below and disagrees with the C++ on any soft-clipped read.
# ---------------------------------------------------------------------------


class MalformedBamError(Exception):
    """The BAM is unreadable past a record htslib will not parse."""


# ---------------------------------------------------------------------------
# MALFORMED_RECORDS
#
# snap-aligner 2.0.5 can emit records whose CIGAR does not consume the whole
# query, which htslib rejects as it parses them:
#
#   [E::bam_read1] CIGAR and query sequence lengths differ for <read>
#   OSError: error -4 while reading file
#
# The cause is alignments at a contig boundary, not the alignment caps. This
# is snap's oldest and most persistent bug family -- its author's own summary
# is "related to things that align right against the end of the contig
# (usually SNAP bugs are)" (snap-user, "CIGAR and query sequence of different
# length"), amplab/snap#121 reports ~1000 such records per 1M reads against a
# small reference, and boundary CIGAR fixes appear in the release notes for
# 2.0.0, 2.0.1, 2.0.2 and again in 2.0.5, the current release ("bugfixes for
# the case where indels would move the start or end of an alignment beyond a
# contig boundary"). Residual cases clearly remain.
#
# That makes it a property of the *assembly*, not of our flags. Measured
# across one ORP library's assemblies: several complete even with the
# per-contig cap disabled, while SRR1789336_80Threads.ORP.fasta fails
# whatever the cap is set to. Loosening -mpc is therefore neither necessary
# nor sufficient to produce it. (An earlier revision of this comment blamed
# snap's secondary-alignment path and told users to set -mpc 1, on the
# strength of a single -mpc 0 run. The cross-assembly comparison falsifies
# that.)
#
# Skipping the record and continuing is safe *for this error*. htslib's
# bam_read1 reads the record's bytes off the BGZF stream in full and only
# then checks the CIGAR against l_qseq (htslib sam.c), so the stream is
# already positioned at the next record when it returns -4. Verified against
# a hand-built BAM carrying malformed records: every subsequent record reads
# back correctly.
#
# But -4 is also what a truncated or corrupt BAM returns, and there the
# stream *is* desynchronised, so continuing would be nonsense. The two are
# indistinguishable from Python -- htslib's reason goes to its own log, not
# into the exception -- so they are told apart by shape instead: snap's bad
# records are isolated, corruption is not. MALFORMED_RUN_LIMIT consecutive
# failures ends the run rather than scoring rubbish.
#
# Skipped alignments are lost coverage, so the count is never swallowed: the
# caller reports it, and the fragments those alignments carried simply go
# unmapped.
# ---------------------------------------------------------------------------

#: Consecutive unparseable records taken as corruption rather than snap's
#: contig-boundary bug. See MALFORMED_RECORDS.
MALFORMED_RUN_LIMIT = 100


@dataclass
class MalformedRecordStats:
    """How many records htslib refused during one pass over a BAM."""

    skipped: int = 0

    def describe(self) -> str:
        return (
            f"{self.skipped} alignment(s) skipped: snap wrote a CIGAR that "
            "does not match the read length, a known snap-aligner bug at "
            "contig boundaries. Those fragments count as unmapped"
        )


def iter_alignments(bam, path="", stats=None):
    """Iterate a BAM, skipping records htslib refuses to parse.

    Args:
        bam: an open :class:`pysam.AlignmentFile`.
        path: the BAM's path, quoted in the error if the file turns out to
            be corrupt rather than merely carrying snap's bad records.
        stats: optional :class:`MalformedRecordStats` to count skips into.
            Callers that want to report them must supply one.

    Raises:
        MalformedBamError: on MALFORMED_RUN_LIMIT consecutive failures, which
            means the BAM is truncated or corrupt.  See MALFORMED_RECORDS.
    """
    iterator = bam.fetch(until_eof=True)
    consecutive = 0
    while True:
        try:
            read = next(iterator)
        except StopIteration:
            return
        except OSError as error:
            consecutive += 1
            if stats is not None:
                stats.skipped += 1
            if consecutive >= MALFORMED_RUN_LIMIT:
                raise MalformedBamError(
                    f"{MALFORMED_RUN_LIMIT} unreadable records in a row"
                    f"{f' in {path}' if path else ''}: {error}\n"
                    "htslib names the reason on the line above, e.g. 'CIGAR "
                    "and query sequence lengths differ'.\n"
                    "Isolated records like that are a known snap-aligner "
                    "bug at contig boundaries and are skipped, but this many "
                    "together means the BAM is truncated or corrupt -- "
                    "usually an aligner run that died partway. Delete it and "
                    "map again."
                ) from error
            continue
        consecutive = 0
        yield read


@dataclass
class ContigMetrics:
    """Accumulator for one reference sequence."""

    name: str
    length: int
    coverage: np.ndarray = field(repr=False, default=None)

    reads_mapped: int = 0
    bases_mapped: int = 0
    p_seq_true_sum: float = 0.0
    fragments_mapped: int = 0
    both_mapped: int = 0
    properpair: int = 0
    bridges: int = 0
    good: int = 0
    bases_uncovered: int = 0
    p_not_segmented: float = 0.0

    # Soft-clip accounting. Kept because it is the single most useful number
    # for interpreting sCcov under snap 2.x, and the BAM it comes from is
    # normally deleted at the end of a run. See SOFT_CLIP_FIX.
    clipped_alignments: int = 0
    clipped_bases: int = 0
    leading_clipped_bases: int = 0

    def __post_init__(self):
        if self.coverage is None:
            self.coverage = np.zeros(max(self.length, 0), dtype=np.int32)

    def add_alignment(self, read: pysam.AlignedSegment) -> None:
        """Add one alignment's coverage. See SOFT_CLIP_FIX."""
        cigar = read.cigartuples
        if cigar is None:
            return
        pos = read.reference_start
        ref_length = self.coverage.size
        clipped = 0
        for index, (op, op_len) in enumerate(cigar):
            if op in _COVERING_OPS:
                end = pos + op_len
                if pos < ref_length:
                    self.coverage[pos : min(end, ref_length)] += 1
                pos = end
            elif op in _SKIPPING_OPS:
                pos += op_len
            elif op == _CIGAR_SOFT_CLIP:
                # Consumes query only -- the cursor must not move. Counted
                # so the run can report how much clipping occurred without
                # anyone having to keep the BAM.
                clipped += op_len
                if index == 0:
                    self.leading_clipped_bases += op_len
            # I consumes query only; H, P consume neither. No cursor move.
        if clipped:
            self.clipped_alignments += 1
            self.clipped_bases += clipped

    def calculate_uncovered_bases(self) -> None:
        self.bases_uncovered = int(np.count_nonzero(self.coverage == 0))

    def set_p_not_segmented(self, nullprior: float = DEFAULT_NULL_PRIOR) -> None:
        states = bin_coverage(self.coverage)
        self.p_not_segmented = prob_not_segmented(states, nullprior=nullprior)

    def p_seq_true(self) -> float:
        if self.reads_mapped > 0:
            return self.p_seq_true_sum / self.reads_mapped
        return 1.0  # bam-read emits a literal 1 for contigs with no reads

    def as_row(self) -> dict:
        return {
            "name": self.name,
            "p_seq_true": self.p_seq_true(),
            "bridges": self.bridges,
            "length": self.length,
            "fragments_mapped": self.fragments_mapped,
            "both_mapped": self.both_mapped,
            "properpair": self.properpair,
            "good": self.good,
            "bases_uncovered": self.bases_uncovered,
            "p_not_segmented": self.p_not_segmented,
        }


def _read_length(read: pysam.AlignedSegment) -> int:
    """Query length, matching BamTools' ``alignment.Length``."""
    if read.query_length:
        return read.query_length
    inferred = read.infer_read_length()
    return inferred if inferred else 0


# ---------------------------------------------------------------------------
# FRAGMENT_ESTIMATOR
#
# estimate_fragment_size() in bam-read.cpp walks the first 10,000 alignments,
# pairs *adjacent* records sharing a name, and accumulates fragment lengths.
# Two quirks are reproduced here verbatim:
#
#   1. The Welford update divides by `count` where the standard recurrence
#      divides by `count + 1` (count is incremented afterwards).  The running
#      mean therefore overshoots and the resulting sd is inflated.
#   2. `mean` is declared `int`, so `mean = mean/(double)count` truncates.
#
# Both feed `realistic_distance = 3*sd + mean`, which gates whether a pair is
# counted as `good` -- and p_good is a direct multiplicand of the contig
# score.  A "corrected" estimator would move every score with nothing to
# validate the new value against, so the published behaviour stands.
# Pass ``faithful=False`` to use a plain two-pass mean and sample sd.
# ---------------------------------------------------------------------------


def estimate_realistic_distance(
    bam_path: str,
    limit: int = 10000,
    faithful: bool = True,
) -> int:
    """Estimate the largest plausible distance between mate start positions.

    Mirrors ``BamRead::estimate_fragment_size``; see FRAGMENT_ESTIMATOR.

    .. important::
       The BAM must be in **read order** (mates adjacent), as aligners emit
       it.  Like the C++, this pairs *adjacent* records sharing a name, so a
       coordinate-sorted BAM separates every mate and yields 0.

    Returns:
        ``3 * sd + mean`` of observed fragment lengths, as an int.  0 when no
        pairs were found.
    """
    fragments: list[int] = []

    name = ""
    prev = ""
    pos1 = pos2 = -1
    len1 = len2 = -1

    count = 0
    running_mean = -1.0
    running_s = 0.0
    total = 0

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        # Tolerant iteration: this pass runs before scoring, and a sample of
        # fragment lengths does not care about a few dropped alignments.
        for read in iter_alignments(bam, bam_path):
            if count >= limit:
                break

            if name != "":
                prev = name
                pos2 = pos1
                len2 = len1

            if read.is_secondary or read.is_supplementary:
                continue

            name = read.query_name or ""
            pos1 = read.reference_start if not read.is_unmapped else -1
            len1 = _read_length(read)

            if prev != name or pos1 < 0 or pos2 < 0:
                continue

            is_reversed = read.is_reverse
            is_mate_reversed = read.mate_is_reverse
            mate_pos = read.next_reference_start

            # Discard pairs pointing away from each other.
            if not is_reversed and is_mate_reversed:
                if pos1 > mate_pos:
                    continue
            elif is_reversed and not is_mate_reversed:
                if mate_pos > pos1:
                    continue

            if pos1 > pos2:
                fragment = pos1 - pos2 + len1
            else:
                fragment = pos2 - pos1 + len2

            fragments.append(fragment)

            if count > 0:
                # Quirk 1: the C++ divides by count, not count + 1.
                new_mean = running_mean + (fragment - running_mean) / count
                running_s += (fragment - running_mean) * (fragment - new_mean)
                running_mean = new_mean
            else:
                running_mean = float(fragment)

            total += fragment
            count += 1

    if count == 0:
        return 0

    if not faithful:
        arr = np.asarray(fragments, dtype=np.float64)
        mean = float(arr.mean())
        sd = float(arr.std(ddof=1)) if arr.size > 1 else 0.0
        return int(3 * sd + mean)

    mean = total // count  # Quirk 2: integer truncation.
    sd = math.sqrt(running_s / (count - 1)) if count > 1 else 0.0
    return int(3 * sd + mean)


def accumulate_metrics(
    references,
    lengths,
    alignments,
    realistic_distance: int,
    nullprior: float = DEFAULT_NULL_PRIOR,
) -> list[ContigMetrics]:
    """Accumulate per-contig metrics from a stream of alignments.

    Takes an iterable rather than a path so that
    :func:`~pytransrate.assign.assign_fragments` can feed alignments straight
    through without staging them in a file -- the round trip salmon 0.8.2's
    ``postSample.bam`` used to force.

    Args:
        references: reference names, in header order.
        lengths: reference lengths, parallel to ``references``.
        alignments: iterable of :class:`pysam.AlignedSegment`, carrying at
            most one alignment per fragment per reference.
        realistic_distance: pairs further apart than this are not ``good``.
        nullprior: prior on "not segmented" passed to the segmenter.

    Returns:
        One :class:`ContigMetrics` per reference, in header order.
    """
    contigs = [
        ContigMetrics(name=name, length=length)
        for name, length in zip(references, lengths)
    ]

    for read in alignments:
        if read.is_unmapped:
            continue

        refid = read.reference_id
        contig = contigs[refid]
        contig.reads_mapped += 1
        contig.bases_mapped += _read_length(read)
        contig.add_alignment(read)

        # Rescaled per-base sequence accuracy from the edit distance.
        # A read with NM == 35 scores 0; the C++ does not clamp, so
        # heavily-mismatched reads contribute negative values.
        if read.has_tag("NM"):
            nm = read.get_tag("NM")
            length = _read_length(read)
            if length > 0:
                scale = (length - 35) / length
                seq_true = (length - nm) / length
                if scale != 1.0:
                    seq_true = (seq_true - scale) * (1 / (1 - scale))
                contig.p_seq_true_sum += seq_true

        is_first = read.is_read1
        is_second = read.is_read2
        mate_mapped = read.is_paired and not read.mate_is_unmapped

        if is_first or (is_second and not mate_mapped):
            contig.fragments_mapped += 1

        # Everything below is counted once per fragment, from read 1.
        if not (is_first and mate_mapped):
            continue

        contig.both_mapped += 1

        if read.is_proper_pair:
            contig.properpair += 1

        if refid != read.next_reference_id:
            contig.bridges += 1
            continue

        mate_pos = read.next_reference_start
        pos = read.reference_start
        if abs(pos - mate_pos) > realistic_distance:
            continue

        # Expect mates on opposite strands, inner-facing.
        if not read.is_reverse and read.mate_is_reverse:
            if pos < mate_pos:
                contig.good += 1
        elif read.is_reverse and not read.mate_is_reverse:
            if mate_pos < pos:
                contig.good += 1

    for contig in contigs:
        contig.calculate_uncovered_bases()
        contig.set_p_not_segmented(nullprior=nullprior)

    return contigs


def compute_bam_metrics(
    bam_path: str,
    nullprior: float = DEFAULT_NULL_PRIOR,
    realistic_distance: int | None = None,
) -> list[ContigMetrics]:
    """Compute per-contig metrics from a BAM on disk.

    Convenience wrapper over :func:`accumulate_metrics` for a BAM that
    already carries one alignment per fragment.  The full pipeline instead
    streams assigned alignments straight in; see
    :func:`~pytransrate.assign.assign_fragments`.
    """
    if realistic_distance is None:
        realistic_distance = estimate_realistic_distance(bam_path)

    stats = MalformedRecordStats()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        metrics = accumulate_metrics(
            bam.references,
            bam.lengths,
            iter_alignments(bam, bam_path, stats),
            realistic_distance=realistic_distance,
            nullprior=nullprior,
        )
    if stats.skipped:
        logger.warning("%s", stats.describe())
    return metrics


def write_metrics_csv(contigs: list[ContigMetrics], path: str) -> None:
    """Write the bam-read-compatible CSV."""
    import csv

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for contig in contigs:
            writer.writerow(contig.as_row())
