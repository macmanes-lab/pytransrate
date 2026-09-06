"""Assigning each fragment to a single transcript.

This replaces salmon 0.8.2's ``--sampleOut``, which wrote a ``postSample.bam``
where every multi-mapping fragment had been sampled down to one alignment.
salmon 2.x accepts the flag but does nothing with it -- it logs "accepted but
not yet implemented and have no effect" -- so the step has moved in-process.

That is a redesign, not a reproduction.  salmon 0.8.2 sampled from the
posterior of its own model (error model, fragment-length distribution, bias
terms); none of that is recoverable from outside.  What follows is a
deterministic maximum-a-posteriori assignment built from what the BAM and
``quant.sf`` actually carry.  It will not reproduce 0.8.2's numbers, and is
not trying to -- see ASSIGNMENT_MODEL.

Because it consumes alignments and yields alignments, the caller can stream
straight into :func:`transrate.bam_metrics.accumulate_metrics` without ever
staging a BAM on disk -- removing a write and a re-read that the
``postSample.bam`` design forced.
"""

from __future__ import annotations

import math
from collections import defaultdict

import pysam

__all__ = [
    "DEFAULT_ERROR_RATE",
    "ORPHAN_EDIT_FRACTION",
    "PRIOR_PSEUDOCOUNT",
    "assign_fragments",
    "group_by_fragment",
    "score_candidates",
]

# ---------------------------------------------------------------------------
# ASSIGNMENT_MODEL
#
# For a fragment f and candidate transcript t we rank by
#
#     log P(t | f)  ~  log prior(t) + sum_a log P(a | t)
#
# prior(t) is salmon's estimated fragment count for t (`NumReads` in
# quant.sf), which is exactly "how many fragments we expect from t" and so is
# the right prior -- TPM alone would ignore transcript length.
#
# P(a | t) comes from the alignment's edit distance:
#
#     log P(a | t) = NM * log(e) + (len - NM) * log(1 - e)
#
# NM rather than AS because snap-aligner emits no AS tag at all (verified on
# 2.0.5: every record carries NM, none carries AS -- salmon warns about this
# too).
#
# Every candidate must explain the *same* observations, or the comparison is
# not a likelihood comparison at all: a candidate carrying one mate would
# otherwise beat a candidate carrying two simply by multiplying fewer terms
# below 1. So a mate missing from a candidate is charged as though it had
# aligned badly -- edit distance ORPHAN_EDIT_FRACTION of its length -- using
# the same formula. That scales with read length automatically, and when a
# mate genuinely aligns nowhere the charge is identical across candidates and
# cancels out of the ranking.
#
# Ties are broken deterministically: higher prior, then lexicographically
# smaller reference name. Never by BAM order, so the result does not depend
# on how the aligner happened to emit records. This matches the direction
# salmon 2.6 took in making deterministic quantification the default.
# ---------------------------------------------------------------------------

#: Per-base sequencing error rate assumed when converting NM to a likelihood.
DEFAULT_ERROR_RATE = 0.01

#: A mate absent from a candidate is scored as if this fraction of its bases
#: mismatched.
#:
#: This sets where a proper pair stops being preferred over a single perfect
#: mate. Measured on 100bp reads at the default error rate, the crossover
#: falls between 7 and 8 mismatches per mate: below that the pair wins, above
#: it the half-placement does. That is the intended shape -- a pair diverging
#: by more than ~7.5% per mate is not evidence the fragment came from that
#: transcript -- but it is a tunable, not a law. Raise it to insist on pairs
#: more strongly.
ORPHAN_EDIT_FRACTION = 0.15

#: Added to salmon's fragment count so transcripts salmon gave zero mass
#: remain reachable rather than being excluded outright.
PRIOR_PSEUDOCOUNT = 1e-3


def group_by_fragment(alignments):
    """Yield ``(name, [alignments])`` for each fragment.

    Relies on the BAM being in **read order**, so a fragment's records --
    both mates and every secondary alignment -- are contiguous.  That is how
    aligners emit BAMs, and was verified on snap-aligner 2.0.5 output.  A
    coordinate-sorted BAM will silently split fragments apart.
    """
    current = None
    batch: list = []
    for read in alignments:
        name = read.query_name
        if name != current:
            if batch:
                yield current, batch
            current = name
            batch = []
        batch.append(read)
    if batch:
        yield current, batch


def _log_likelihood(nm: float, length: int, error_rate: float) -> float:
    """Log P(alignment | transcript) for `nm` mismatches over `length` bases."""
    if length <= 0:
        return 0.0
    nm = max(0.0, min(float(nm), float(length)))
    return nm * math.log(error_rate) + (length - nm) * math.log1p(-error_rate)


def _read_length(read) -> int:
    return read.query_length or read.infer_read_length() or 0


def _log_alignment_likelihood(read, error_rate: float) -> float:
    """Log P(alignment | transcript), from the edit distance."""
    length = _read_length(read)
    try:
        nm = int(read.get_tag("NM"))
    except KeyError:
        nm = 0
    return _log_likelihood(nm, length, error_rate)


def score_candidates(
    batch,
    references,
    priors,
    error_rate: float = DEFAULT_ERROR_RATE,
    orphan_edit_fraction: float = ORPHAN_EDIT_FRACTION,
):
    """Score each transcript this fragment could have come from.

    Args:
        batch: the fragment's alignments.
        references: reference names by id.
        priors: reference name -> expected fragment count.
        error_rate: per-base error rate for the likelihood.
        orphan_edit_fraction: mismatch fraction charged for a mate missing
            from a candidate. See ASSIGNMENT_MODEL.

    Returns:
        ``{reference_id: (score, prior, name)}``.
    """
    by_ref: dict[int, list] = defaultdict(list)
    mates_present = set()
    typical_length = 0
    for read in batch:
        if read.is_unmapped:
            continue
        by_ref[read.reference_id].append(read)
        mates_present.add(bool(read.is_read2))
        typical_length = max(typical_length, _read_length(read))

    if not mates_present:
        return {}

    # Charge a missing mate as a badly-aligned one, so every candidate is
    # scored over the same set of mates.
    orphan_cost = _log_likelihood(
        orphan_edit_fraction * typical_length, typical_length, error_rate
    )

    scored = {}
    for ref_id, reads in by_ref.items():
        name = references[ref_id]
        prior = priors.get(name, 0.0) + PRIOR_PSEUDOCOUNT

        explained = {bool(r.is_read2) for r in reads}
        score = math.log(prior)
        for read in reads:
            score += _log_alignment_likelihood(read, error_rate)
        score += orphan_cost * len(mates_present - explained)

        scored[ref_id] = (score, prior, name)
    return scored


def assign_fragments(
    alignments,
    references,
    expression=None,
    error_rate: float = DEFAULT_ERROR_RATE,
    orphan_edit_fraction: float = ORPHAN_EDIT_FRACTION,
    clear_secondary: bool = True,
):
    """Yield one transcript's worth of alignments per fragment.

    Args:
        alignments: read-ordered alignments (see :func:`group_by_fragment`).
        references: reference names by id, e.g. ``bam.references``.
        expression: salmon output from
            :func:`transrate.quantify.load_expression`.  When omitted every
            transcript gets an equal prior and assignment falls back to
            alignment quality alone.
        error_rate: per-base error rate for the likelihood.
        orphan_edit_fraction: mismatch fraction charged per mate missing
            from a candidate.
        clear_secondary: unset the secondary flag on emitted records, since
            the surviving alignment is now the fragment's only placement.

    Yields:
        :class:`pysam.AlignedSegment`, in input order within each fragment.
    """
    priors = {}
    if expression:
        for name, values in expression.items():
            priors[name] = float(values.get("eff_count", 0.0))

    for _name, batch in group_by_fragment(alignments):
        scored = score_candidates(
            batch, references, priors, error_rate, orphan_edit_fraction
        )
        if not scored:
            continue

        # Deterministic: best score, then higher prior, then name.
        best_id = min(
            scored,
            key=lambda ref_id: (
                -scored[ref_id][0],
                -scored[ref_id][1],
                scored[ref_id][2],
            ),
        )

        for read in batch:
            if read.is_unmapped or read.reference_id != best_id:
                continue
            if clear_secondary and read.is_secondary:
                read.flag = read.flag & ~pysam.FSECONDARY
            yield read
