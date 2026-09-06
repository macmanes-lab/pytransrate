"""Sequence statistics: base composition and longest ORF.

Port of ``ext/transrate/transrate.c``, the C extension the Ruby transrate
compiled and mixed into ``Transrate::Contig``.  ``kmer_count`` from that file
is not ported -- nothing outside the Ruby test suite ever called it.

Both functions reproduce the C semantics exactly, including its case
handling, which differs between the two:

* :func:`base_composition` folds case (``if (base > 90) base -= 32``), so
  lowercase sequence is counted normally.
* :func:`longest_orf` compares against uppercase literals only, so a
  soft-masked (lowercase) contig yields an ORF length of 0.  See
  ORF_CASE_SENSITIVITY.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "BASES",
    "base_composition",
    "dibase_composition",
    "longest_orf",
]

#: Base order used by the C extension's index mapping (A=0, C=1, G=2, T=3, N=4).
BASES = ("a", "c", "g", "t", "n")

_A, _C, _G, _T = (ord(c) for c in "ACGT")

# Maps every byte to its index in BASES, after case folding.
_BASE_INDEX = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate("ACGT"):
    _BASE_INDEX[ord(_ch)] = _i
    _BASE_INDEX[ord(_ch.lower())] = _i


def _as_bytes(seq) -> np.ndarray:
    if isinstance(seq, str):
        seq = seq.encode("ascii", errors="replace")
    return np.frombuffer(bytes(seq), dtype=np.uint8)


def base_composition(seq) -> dict[str, int]:
    """Count A/C/G/T/N, case-insensitively.

    Returns:
        Counts keyed by ``BASES``; anything that is not ACGT counts as ``n``.
    """
    arr = _as_bytes(seq)
    if arr.size == 0:
        return {b: 0 for b in BASES}
    idx = _BASE_INDEX[arr]
    counts = np.bincount(idx, minlength=5)
    return {b: int(counts[i]) for i, b in enumerate(BASES)}


def dibase_composition(seq) -> dict[str, int]:
    """Count the 25 ordered dinucleotides, case-insensitively.

    Mirrors the C indexing ``prev * 5 + curr`` over overlapping pairs.
    """
    arr = _as_bytes(seq)
    keys = [f"{a}{b}" for a in BASES for b in BASES]
    if arr.size < 2:
        return {k: 0 for k in keys}
    idx = _BASE_INDEX[arr].astype(np.int64)
    pair = idx[:-1] * 5 + idx[1:]
    counts = np.bincount(pair, minlength=25)
    return {k: int(counts[i]) for i, k in enumerate(keys)}


# ---------------------------------------------------------------------------
# ORF_CASE_SENSITIVITY
#
# method_longest_orf in transrate.c compares raw bytes against 'A', 'T', 'G'
# and never folds case, unlike method_composition immediately above it in the
# same file.
#
# The consequence is the opposite of the intuitive one.  Lowercase input does
# not yield an ORF length of 0 -- it yields the *maximum* possible length.  No
# stop codon is ever recognised, so the frame counter never resets and the
# whole contig reads as a single open frame.  Measured against the C for a
# random 3 kb sequence: 90 codons uppercase, 1000 codons lowercase.
#
# So a soft-masked assembly silently inflates n_with_orf and
# mean_orf_percent rather than obviously zeroing them.  Preserved anyway:
# orf_length feeds only those two assembly statistics and the contigs.csv
# column of the same name -- never the contig score -- so nothing in ORP's
# contig selection depends on it, and matching the C keeps the assembly-level
# statistics comparable with published runs.  Upper-case the sequence before
# calling if you need ORFs from soft-masked input.
# ---------------------------------------------------------------------------


def _classify_codons(arr: np.ndarray):
    """Boolean masks for forward start/stop codons at each offset."""
    a, b, c = arr[:-2], arr[1:-1], arr[2:]
    is_start = (a == _A) & (b == _T) & (c == _G)
    is_stop = (a == _T) & (
        ((b == _A) & (c == _G))  # TAG
        | ((b == _A) & (c == _A))  # TAA
        | ((b == _G) & (c == _A))  # TGA
    )
    return is_start, is_stop


def _classify_codons_reverse(arr: np.ndarray):
    """Masks for the reverse strand, read on the forward sequence.

    The C walks i downward testing ``str[i], str[i-1], str[i-2]``; at triplet
    offset ``j = i - 2`` that is CAT for a start and CTA/TTA/TCA for a stop --
    the reverse complements of ATG and TAG/TAA/TGA.
    """
    a, b, c = arr[:-2], arr[1:-1], arr[2:]
    is_start = (a == _C) & (b == _A) & (c == _T)
    is_stop = (c == _A) & (
        ((a == _C) & (b == _T))  # CTA
        | ((a == _T) & (b == _T))  # TTA
        | ((a == _T) & (b == _C))  # TCA
    )
    return is_start, is_stop


def _longest_in_frame(is_start: np.ndarray, is_stop: np.ndarray) -> int:
    """Run the C's per-frame state machine over one frame's codons.

    The C keeps ``len[i%3]``, incrementing on every codon while non-negative,
    zeroing to -1 at a stop, and restarting at 1 on the next start.  So the
    frame's best run is: codons before the first stop, then for each later
    segment, the distance from its first start to the segment end.
    """
    m = int(is_start.size)
    if m == 0:
        return 0

    stop_idx = np.flatnonzero(is_stop)
    if stop_idx.size == 0:
        # len starts at 0 and never resets, so the whole frame counts.
        return m

    # From the sequence start up to the first stop, every codon counts --
    # len begins at 0, which is >= 0, so no start codon is required.
    best = int(stop_idx[0])

    start_idx = np.flatnonzero(is_start)
    for k in range(stop_idx.size):
        seg_begin = int(stop_idx[k]) + 1
        seg_end = int(stop_idx[k + 1]) if k + 1 < stop_idx.size else m
        # Only the first start in the segment matters; later ones merely
        # increment an already-open reading frame.
        p = int(np.searchsorted(start_idx, seg_begin, side="left"))
        if p < start_idx.size and int(start_idx[p]) < seg_end:
            best = max(best, seg_end - int(start_idx[p]))
    return best


def longest_orf(seq) -> int:
    """Longest open reading frame across all six frames, in codons.

    An ORF runs from either the sequence start or a start codon, to either a
    stop codon or the sequence end.  Uppercase only; see
    ORF_CASE_SENSITIVITY.

    Returns:
        Length in **codons**, matching the C.  ``basic_stats`` multiplies by
        3 to reach bases.
    """
    arr = _as_bytes(seq)
    n = int(arr.size)
    if n < 3:
        return 0

    longest = 0

    fwd_start, fwd_stop = _classify_codons(arr)
    for frame in range(3):
        longest = max(
            longest, _longest_in_frame(fwd_start[frame::3], fwd_stop[frame::3])
        )

    rev_start, rev_stop = _classify_codons_reverse(arr)
    offsets = np.arange(rev_start.size)
    for frame in range(3):
        # The C indexes len[i%3] with i = j + 2, and walks j downward.
        sel = offsets[(offsets + 2) % 3 == frame][::-1]
        longest = max(longest, _longest_in_frame(rev_start[sel], rev_stop[sel]))

    return int(longest)
