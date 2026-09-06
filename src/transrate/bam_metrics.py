"""Per-contig metrics from a BAM alignment.

Replaces the external ``bam-read`` binary (Blahah/transrate-tools), which the
Ruby transrate shelled out to.  Produces the same ten columns:

    name, p_seq_true, bridges, length, fragments_mapped, both_mapped,
    properpair, good, bases_uncovered, p_not_segmented

One deliberate behavioural change from the C++ is documented at SOFT_CLIP_FIX.
Two quirks are preserved on purpose and documented at FRAGMENT_ESTIMATOR and
in :mod:`transrate.segmenter` (BINNING_QUIRK).

The input is expected to carry exactly one alignment per fragment -- the role
``postSample.bam`` played under salmon 0.8.2's ``--sampleOut``.  See
:mod:`transrate.assign`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pysam

from transrate.segmenter import DEFAULT_NULL_PRIOR, bin_coverage, prob_not_segmented

__all__ = [
    "CSV_COLUMNS",
    "ContigMetrics",
    "estimate_realistic_distance",
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
        for op, op_len in cigar:
            if op in _COVERING_OPS:
                end = pos + op_len
                if pos < ref_length:
                    self.coverage[pos : min(end, ref_length)] += 1
                pos = end
            elif op in _SKIPPING_OPS:
                pos += op_len
            # I, S consume query only; H, P consume neither. No cursor move.

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
        for read in bam.fetch(until_eof=True):
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


def compute_bam_metrics(
    bam_path: str,
    nullprior: float = DEFAULT_NULL_PRIOR,
    realistic_distance: int | None = None,
) -> list[ContigMetrics]:
    """Compute per-contig metrics, one entry per reference in the BAM header.

    Args:
        bam_path: BAM carrying one alignment per fragment.
        nullprior: prior on "not segmented" passed to the segmenter.
        realistic_distance: override the fragment-size estimate; computed from
            the BAM when omitted.

    Returns:
        One :class:`ContigMetrics` per reference, in header order.
    """
    if realistic_distance is None:
        realistic_distance = estimate_realistic_distance(bam_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contigs = [
            ContigMetrics(name=name, length=length)
            for name, length in zip(bam.references, bam.lengths)
        ]

        for read in bam.fetch(until_eof=True):
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


def write_metrics_csv(contigs: list[ContigMetrics], path: str) -> None:
    """Write the bam-read-compatible CSV."""
    import csv

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for contig in contigs:
            writer.writerow(contig.as_row())
