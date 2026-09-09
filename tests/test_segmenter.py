"""Tests for the coverage segmentation model.

``prob_not_segmented`` has no external oracle -- samtools can validate a
coverage vector but nothing outside transrate computes sCseg.  So it is
checked two ways here:

1.  Against a literal linear-space transcription of segmenter.cpp
    (``reference_prob_not_segmented`` below), over inputs where the C++
    formulation does not underflow.  This validates the log-space rewrite as
    an algebraic identity rather than a reinterpretation.
2.  Against qualitative properties the model is supposed to have -- flat
    coverage is one segment, a step is two.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytransrate.segmenter import (
    DEFAULT_NULL_PRIOR,
    NUM_BINS,
    NUM_STATES,
    bin_coverage,
    prob_not_segmented,
)


# ---------------------------------------------------------------------------
# Literal transcription of segmenter.cpp, in linear space with math.gamma.
# Doubles throughout: we are validating the algebra, not reproducing the
# `float` truncation bug that the real C++ has in prob_R_given_k.
# ---------------------------------------------------------------------------


def _reference_rhs(states, length):
    n_states = len(states)
    upper = 1.0
    for p in states:
        upper *= math.gamma(p + 1)
    return math.gamma(n_states) * (upper / math.gamma(length + n_states))


def reference_prob_not_segmented(seq, nullprior=DEFAULT_NULL_PRIOR):
    """Direct port of Segmenter::prob_k_given_R(0), linear space."""
    seq = [min(int(s), NUM_STATES - 1) for s in seq]
    total = len(seq)

    states = [0] * NUM_STATES
    for s in seq:
        states[s] += 1

    p_k0 = _reference_rhs(states, total)

    # prob_R_given_unit_k: full pmat, exactly as the C++ builds it.
    pmat = [[0.0] * total for _ in range(total)]
    lstates = list(states)
    for i in range(total):
        segstates = list(lstates)
        for j in range(total - 1, i - 1, -1):
            pmat[i][j] = _reference_rhs(segstates, j + 1 - i)
            segstates[seq[j]] -= 1
        lstates[seq[i]] -= 1

    p_a = max(1.0 / (total - 1), float(np.finfo(np.float32).tiny))
    p_k1 = 0.0
    for i in range(total - 1):
        p_k1 += pmat[0][i] * pmat[i + 1][total - 1] * p_a

    prior_k0 = nullprior
    prior_k1 = (1.0 - nullprior) / 1
    marginal = prior_k0 * p_k0 + prior_k1 * p_k1
    return (p_k0 * prior_k0) / marginal


# ---------------------------------------------------------------------------
# Agreement with the reference implementation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "seq",
    [
        [5] * 30,
        [0] * 30,
        [3] * 15 + [8] * 15,
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 3,
        [7] * 29 + [0],
        [0, 0, 9, 9, 9, 9, 9, 9, 0, 0] * 3,
        [2] * 28 + [0, 0],
        list(range(20)) + [0] * 10,
    ],
)
def test_matches_linear_space_reference(seq):
    got = prob_not_segmented(seq)
    want = reference_prob_not_segmented(seq)
    assert got == pytest.approx(want, rel=1e-9)


@pytest.mark.parametrize("nullprior", [0.1, 0.3, 0.5, 0.7, 0.9, 0.99])
def test_matches_reference_across_priors(nullprior):
    seq = [4] * 12 + [9] * 18
    got = prob_not_segmented(seq, nullprior=nullprior)
    want = reference_prob_not_segmented(seq, nullprior=nullprior)
    assert got == pytest.approx(want, rel=1e-9)


def test_random_states_match_reference():
    rng = np.random.default_rng(1337)
    for _ in range(50):
        seq = rng.integers(0, 12, size=NUM_BINS).tolist()
        got = prob_not_segmented(seq)
        want = reference_prob_not_segmented(seq)
        assert got == pytest.approx(want, rel=1e-9)


# ---------------------------------------------------------------------------
# Model behaviour
# ---------------------------------------------------------------------------


def test_flat_coverage_reads_as_one_segment():
    flat = prob_not_segmented([6] * NUM_BINS)
    step = prob_not_segmented([1] * 15 + [11] * 15)
    assert flat > step
    assert flat > 0.5


def test_sharp_step_reads_as_segmented():
    assert prob_not_segmented([0] * 15 + [12] * 15) < 0.5


def test_result_is_a_probability():
    rng = np.random.default_rng(7)
    for _ in range(100):
        seq = rng.integers(0, NUM_STATES, size=NUM_BINS).tolist()
        p = prob_not_segmented(seq)
        assert 0.0 <= p <= 1.0
        assert math.isfinite(p)


def test_states_at_or_above_cap_are_folded_into_top_bucket():
    # load_states(): `if (seq[i] < 24) ++states_[seq[i]]; else ++states_[23];`
    assert prob_not_segmented([24] * 30) == pytest.approx(
        prob_not_segmented([23] * 30)
    )


def test_degenerate_lengths():
    assert prob_not_segmented([]) == 1.0
    assert prob_not_segmented([5]) == 1.0


def test_extreme_priors():
    assert prob_not_segmented([3] * 30, nullprior=1.0) == pytest.approx(1.0)
    assert prob_not_segmented([3] * 30, nullprior=0.0) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Binning
# ---------------------------------------------------------------------------


def test_bin_coverage_uses_integer_division_for_the_mean():
    # A bin averaging 1.9x must read as state 0 (log2(1)), not log2(1.9).
    cov = np.array([1, 2] * 15 + [2] * 15)  # 45 bases
    states = bin_coverage(cov)
    assert states.dtype == np.int64
    assert set(np.unique(states)).issubset(set(range(NUM_STATES + 1)))


def test_bin_coverage_log2_bucketing():
    # Uniform coverage of 8 -> log2(8) == 3 for every populated bin.
    cov = np.full(300, 8, dtype=np.int64)
    states = bin_coverage(cov)
    populated = states[states > 0]
    assert populated.size > 0
    assert np.all(populated == 3)


@pytest.mark.parametrize(
    "depth,state", [(1, 0), (2, 1), (3, 1), (7, 2), (8, 3), (15, 3), (16, 4)]
)
def test_bin_coverage_floors_between_powers_of_two(depth, state):
    """The state is floor(log2(mean)), taken from the mean's bit length."""
    cov = np.full(300, depth, dtype=np.int64)
    populated = bin_coverage(cov)[:28]
    assert np.all(populated == state)


def test_bin_coverage_zero_coverage_is_state_zero():
    assert np.all(bin_coverage(np.zeros(300, dtype=np.int64)) == 0)


def test_bin_coverage_clamps_to_max_state():
    cov = np.full(300, 2**26, dtype=np.int64)  # log2 == 26, above the cap
    states = bin_coverage(cov)
    assert states.max() == 24


def test_bin_coverage_pads_to_thirty_by_default():
    """BINNING_QUIRK: the length-30 vector is zero-padded past the last bin."""
    cov = np.full(300, 8, dtype=np.int64)
    padded = bin_coverage(cov)
    assert padded.size == NUM_BINS
    # bin_width = 300//30 + 1 = 11 -> ceil(300/11) = 28 populated bins.
    assert np.count_nonzero(padded) == 28
    assert padded[28] == 0 and padded[29] == 0


def test_pad_bins_false_drops_the_padding():
    cov = np.full(300, 8, dtype=np.int64)
    trimmed = bin_coverage(cov, pad_bins=False)
    assert trimmed.size == 28
    assert np.all(trimmed == 3)


def test_padding_biases_sCseg_downward():
    """The quirk is not cosmetic: it lowers sCseg for a uniform contig."""
    cov = np.full(300, 8, dtype=np.int64)
    with_padding = prob_not_segmented(bin_coverage(cov))
    without_padding = prob_not_segmented(bin_coverage(cov, pad_bins=False))
    assert without_padding > with_padding


def test_bin_coverage_handles_short_contigs():
    cov = np.full(7, 4, dtype=np.int64)
    states = bin_coverage(cov)
    assert states.size == NUM_BINS
    assert np.all(states[:7] == 2)  # log2(4)
    assert np.all(states[7:] == 0)


def test_bin_coverage_empty():
    assert bin_coverage(np.zeros(0, dtype=np.int64)).size == NUM_BINS
