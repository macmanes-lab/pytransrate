"""Bayesian single-changepoint model over binned contig coverage.

This is a port of ``src/segmenter.cpp`` and the ``setPNotSegmented`` half of
``src/pileup.cpp`` from Blahah/transrate-tools, which the Ruby transrate
invoked as the external ``bam-read`` binary.

The model asks whether a contig's coverage profile is better explained by one
segment (k=0) or by two segments split at a changepoint (k=1).  ``sCseg`` in
the transrate score is the posterior probability of k=0.

Per-base coverage is binned, each bin is reduced to an integer state
(``floor(log2(mean coverage))``, clamped), and the state vector is scored
under a Dirichlet-multinomial:

    P(R | segment) = Gamma(S) * prod_i Gamma(c_i + 1) / Gamma(n + S)

for S states, per-state counts c_i, and segment length n.  For k=1 this is
marginalised over every changepoint under a flat prior.

Two deliberate departures from the C++ are documented at BINNING_QUIRK and in
``prob_not_segmented``; both are called out where they occur.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln, logsumexp

__all__ = [
    "NUM_STATES",
    "NUM_BINS",
    "MAX_STATE",
    "DEFAULT_NULL_PRIOR",
    "bin_coverage",
    "prob_not_segmented",
]

#: Number of Dirichlet-multinomial categories (``states_`` in segmenter.cpp).
NUM_STATES = 24

#: Number of coverage bins (``vector<int> states(30)`` in pileup.cpp).
NUM_BINS = 30

#: Ceiling applied to a bin's log2 state before bucketing (pileup.cpp).
MAX_STATE = 24

#: Prior on k=0, i.e. "not segmented" (``nullprior`` in segmenter.cpp).
DEFAULT_NULL_PRIOR = 0.7

#: Only k in {0, 1} is considered (``const int maxk = 1``).
_MAX_K = 1

# ---------------------------------------------------------------------------
# BINNING_QUIRK
#
# pileup.cpp declares `vector<int> states(30)` -- zero-initialised, length 30
# unconditionally -- then fills only as many entries as there are bins:
#
#     int bin_width = (ref_length / 30) + 1;
#
# Because bin_width is rounded *up*, the number of bins actually written is
# ceil(ref_length / bin_width), which is strictly less than 30 for every
# ref_length that is not a multiple of 30.  The unwritten tail stays 0, and
# the full length-30 vector is what gets handed to the Segmenter.
#
# State 0 means "uncovered", so those padding entries look to the model like
# an uncovered stretch at the 3' end of every contig.  That biases sCseg
# downward roughly uniformly across an assembly.
#
# We reproduce it.  It is the behaviour of the published, cited method, and
# unlike the soft-clip handling in pileup.py there is no external oracle that
# would tell us a change is an improvement rather than just a difference.
# Set ``pad_bins=False`` to drop the padding and score only real bins.
# ---------------------------------------------------------------------------


def bin_coverage(
    coverage,
    n_bins: int = NUM_BINS,
    max_state: int = MAX_STATE,
    pad_bins: bool = True,
) -> np.ndarray:
    """Reduce a per-base coverage vector to integer states, one per bin.

    Mirrors ``TransratePileup::setPNotSegmented``, including its use of
    integer division for the per-bin mean (so a bin averaging 1.9x reads as
    1x, not 2x).

    Args:
        coverage: per-base coverage counts.
        n_bins: number of bins; also the padded output length.
        max_state: ceiling applied before bucketing.
        pad_bins: keep the zero padding described at BINNING_QUIRK.  True
            reproduces transrate-tools; False scores only populated bins.

    Returns:
        Integer state vector, length ``n_bins`` when ``pad_bins`` is True.
    """
    coverage = np.asarray(coverage, dtype=np.int64)
    ref_length = int(coverage.size)

    states = np.zeros(n_bins, dtype=np.int64)
    if ref_length == 0:
        return states if pad_bins else states[:0]

    # (ref_length / 30) + 1 in C++ integer arithmetic. Rounding up means the
    # bin count never exceeds n_bins, so no clamping is needed below.
    bin_width = ref_length // n_bins + 1

    # Vectorised equivalent of the C++ per-base accumulation loop: sum each
    # bin, then take the integer mean exactly as `total / counter` did.
    starts = np.arange(0, ref_length, bin_width)
    sums = np.add.reduceat(coverage, starts)
    counts = np.diff(np.append(starts, ref_length))
    means = sums // counts  # integer division, as in C++

    # log2(0) is -inf; the C++ max(0.0, ...) pins it to 0. Truncation toward
    # zero matches the C++ (int) cast, and equals floor for means >= 1.
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log2(np.maximum(means, 1))
    binned = np.minimum(max_state, logs.astype(np.int64))
    binned[means <= 0] = 0

    n_filled = binned.size
    states[:n_filled] = binned
    return states if pad_bins else states[:n_filled]


#: gammaln lookup for the denominator, covering the usual bin/state counts.
#: Sized generously; anything beyond falls back to a direct call.
_GAMMALN_TABLE = gammaln(np.arange(4 * (NUM_BINS + NUM_STATES)))


def _log_denominator(length: int, n_states: int) -> float:
    """gammaln(length + n_states), from the table when it reaches."""
    index = length + n_states
    if 0 <= index < _GAMMALN_TABLE.size:
        return float(_GAMMALN_TABLE[index])
    return float(gammaln(index))


def _log_segment_likelihood(counts: np.ndarray, length: int) -> float:
    """Log of ``prob_R_given_k_rhs``: the Dirichlet-multinomial term.

    The C++ evaluates this with ``tgamma`` in linear space and stores the
    result through a ``float`` temporary (``float result = pRk_[k];`` in
    ``prob_R_given_k``).  With 24 categories the numerator and denominator
    reach ~1e22 and ~1e69 respectively, so sparse state vectors produce
    values below FLT_MIN that flush to zero -- and a zero marginal makes the
    posterior 0/0.  Working in log space removes that failure mode entirely;
    it changes the answer only where the C++ had already lost the value.
    """
    n_states = counts.size
    return float(
        gammaln(n_states) + gammaln(counts + 1.0).sum() - gammaln(length + n_states)
    )


def prob_not_segmented(
    states,
    nullprior: float = DEFAULT_NULL_PRIOR,
    n_states: int = NUM_STATES,
) -> float:
    """Posterior probability of k=0 (a single segment) given the states.

    This is ``sCseg`` / ``p_not_segmented``.  Equivalent to
    ``Segmenter::prob_k_given_R(0)`` but evaluated in log space.

    Args:
        states: integer state vector from :func:`bin_coverage`.  Values at or
            above ``n_states`` are folded into the top bucket, matching
            ``load_states``.
        nullprior: prior on k=0.
        n_states: number of Dirichlet-multinomial categories.

    Returns:
        Probability in [0, 1].
    """
    seq = np.asarray(states, dtype=np.int64)
    seq = np.clip(seq, 0, n_states - 1)
    total = int(seq.size)

    if total == 0:
        return 1.0

    counts = np.bincount(seq, minlength=n_states).astype(np.float64)
    log_p_k0 = _log_segment_likelihood(counts, total)

    if total < 2:
        # No changepoint is expressible, so k=1 has no support.
        log_p_k1 = -np.inf
    else:
        # The C++ builds a full total x total matrix but reads only
        # pmat[0][i] (prefix) and pmat[i+1][total-1] (suffix).  Running
        # counts give the same values in O(total * n_states).
        # sum(gammaln(c_i + 1)) changes by exactly log(c + 1) when one count
        # goes c -> c + 1, because gammaln(c + 2) - gammaln(c + 1) = log(c + 1).
        # So the numerator is carried as a running scalar instead of being
        # recomputed over all n_states categories at every step.
        log_gamma_states = float(gammaln(n_states))

        prefix = np.empty(total - 1, dtype=np.float64)
        running = np.zeros(n_states, dtype=np.float64)
        numerator = 0.0
        for i in range(total - 1):
            state = seq[i]
            running[state] += 1.0
            numerator += math.log(running[state])
            prefix[i] = (
                log_gamma_states + numerator - _log_denominator(i + 1, n_states)
            )

        suffix = np.empty(total - 1, dtype=np.float64)
        running = np.zeros(n_states, dtype=np.float64)
        numerator = 0.0
        for i in range(total - 1, 0, -1):
            state = seq[i]
            running[state] += 1.0
            numerator += math.log(running[state])
            suffix[i - 1] = (
                log_gamma_states + numerator - _log_denominator(total - i, n_states)
            )

        # Flat prior over changepoints, floored at FLT_MIN as in the C++.
        p_a = max(1.0 / (total - 1), float(np.finfo(np.float32).tiny))
        log_p_k1 = float(logsumexp(prefix + suffix + math.log(p_a)))

    # Posterior over k in {0, 1}.
    with np.errstate(divide="ignore"):
        log_prior_k0 = math.log(nullprior) if nullprior > 0.0 else -np.inf
        remaining = (1.0 - nullprior) / _MAX_K
        log_prior_k1 = math.log(remaining) if remaining > 0.0 else -np.inf

    log_num_k0 = log_prior_k0 + log_p_k0
    log_num_k1 = log_prior_k1 + log_p_k1

    if math.isinf(log_num_k0) and math.isinf(log_num_k1):
        # Both hypotheses are impossible; fall back to the prior.
        return float(nullprior)

    log_marginal = float(logsumexp([log_num_k0, log_num_k1]))
    return float(math.exp(log_num_k0 - log_marginal))
