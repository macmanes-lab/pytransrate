"""Tests for the C-extension replacements.

The strong check here is differential: ``tests/oracle/orf_oracle.c`` is the
body of ``method_composition`` and ``method_longest_orf`` lifted verbatim out
of ``ext/transrate/transrate.c`` with the Ruby wrappers stripped.  It is
compiled on demand and both implementations are run over the same random
sequences.  If no C compiler is available those tests skip; the hand-written
cases below still run.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pytransrate.sequence import (
    BASES,
    base_composition,
    dibase_composition,
    longest_orf,
)

ORACLE_SRC = Path(__file__).parent / "oracle" / "orf_oracle.c"


@pytest.fixture(scope="session")
def oracle(tmp_path_factory):
    """Compile the extracted C and return a callable over sequences."""
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        pytest.skip("no C compiler available to build the oracle")

    binary = tmp_path_factory.mktemp("oracle") / "orf_oracle"
    result = subprocess.run(
        [cc, "-O2", "-o", str(binary), str(ORACLE_SRC)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.skip(f"oracle failed to build: {result.stderr[:400]}")

    def run(sequences):
        payload = "\n".join(sequences) + "\n"
        proc = subprocess.run(
            [str(binary)], input=payload, capture_output=True, text=True, check=True
        )
        rows = []
        for line in proc.stdout.splitlines():
            a, c, g, t, n, orf = (int(x) for x in line.split())
            rows.append(
                ({"a": a, "c": c, "g": g, "t": t, "n": n}, orf)
            )
        return rows

    return run


def _random_sequences(rng, count, min_len, max_len, alphabet="ACGT"):
    seqs = []
    for _ in range(count):
        n = int(rng.integers(min_len, max_len + 1))
        seqs.append("".join(rng.choice(list(alphabet), size=n)))
    return seqs


# ---------------------------------------------------------------------------
# Differential against the original C
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "alphabet",
    ["ACGT", "ACGTN", "ACGTNacgtn", "ACGTRYKMSWBDHVN"],
)
def test_matches_c_extension(oracle, alphabet):
    rng = np.random.default_rng(20260906)
    seqs = _random_sequences(rng, 200, 1, 400, alphabet)
    expected = oracle(seqs)
    assert len(expected) == len(seqs)
    for seq, (want_comp, want_orf) in zip(seqs, expected):
        assert base_composition(seq) == want_comp, seq
        assert longest_orf(seq) == want_orf, seq


def test_matches_c_extension_on_orf_rich_sequences(oracle):
    """Bias the input toward start/stop codons so the state machine is
    exercised hard rather than incidentally."""
    rng = np.random.default_rng(4242)
    fragments = ["ATG", "TAA", "TAG", "TGA", "CAT", "CTA", "TTA", "TCA", "AAA", "GGC"]
    seqs = []
    for _ in range(300):
        n = int(rng.integers(1, 60))
        seqs.append("".join(rng.choice(fragments, size=n)))
    for seq, (want_comp, want_orf) in zip(seqs, oracle(seqs)):
        assert base_composition(seq) == want_comp, seq
        assert longest_orf(seq) == want_orf, seq


def test_matches_c_extension_on_short_sequences(oracle):
    seqs = []
    for n in range(1, 13):
        seqs.extend(
            "".join(p)
            for p in np.array(
                np.meshgrid(*[list("ATG")] * min(n, 4))
            ).T.reshape(-1, min(n, 4))
        )
    seqs = sorted(set(seqs))
    for seq, (want_comp, want_orf) in zip(seqs, oracle(seqs)):
        assert base_composition(seq) == want_comp, seq
        assert longest_orf(seq) == want_orf, seq


def test_matches_c_extension_on_long_sequences(oracle):
    rng = np.random.default_rng(99)
    seqs = _random_sequences(rng, 12, 5000, 20000, "ACGT")
    for seq, (want_comp, want_orf) in zip(seqs, oracle(seqs)):
        assert base_composition(seq) == want_comp, seq
        assert longest_orf(seq) == want_orf, seq


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def test_composition_counts():
    assert base_composition("AACCCGGGGTTTTT") == {
        "a": 2, "c": 3, "g": 4, "t": 5, "n": 0
    }


def test_composition_folds_case():
    assert base_composition("acgt") == base_composition("ACGT")


def test_composition_treats_ambiguity_codes_as_n():
    comp = base_composition("ACGTNRYKM")
    assert comp["n"] == 5


def test_composition_of_empty_sequence():
    assert base_composition("") == {b: 0 for b in BASES}


def test_dibase_composition_counts_overlapping_pairs():
    dibases = dibase_composition("ACGT")
    assert dibases["ac"] == 1
    assert dibases["cg"] == 1
    assert dibases["gt"] == 1
    assert sum(dibases.values()) == 3


def test_dibase_composition_is_case_insensitive():
    assert dibase_composition("acgt") == dibase_composition("ACGT")


def test_dibase_composition_too_short():
    assert sum(dibase_composition("A").values()) == 0


# ---------------------------------------------------------------------------
# ORF
# ---------------------------------------------------------------------------


def test_orf_counts_codons_not_bases():
    # ATG AAA CCC TAG -> 4 codons before the stop is consumed.
    assert longest_orf("ATGAAACCCTAG") == 4


def test_orf_runs_to_sequence_end_without_a_stop():
    assert longest_orf("ATG" * 10) == 10


def test_orf_restarts_at_a_start_codon_after_a_stop():
    # A stop kills the frame; the next ATG revives it.
    seq = "TAA" + "CCC" * 3 + "ATG" + "AAA" * 5
    assert longest_orf(seq) >= 6


def test_orf_of_short_sequence_is_zero():
    assert longest_orf("") == 0
    assert longest_orf("AT") == 0


def test_orf_finds_reverse_strand():
    # CAT is the reverse complement of ATG; CTA/TTA/TCA are reverse stops.
    forward_only = "CAT" * 8
    assert longest_orf(forward_only) > 0


def test_lowercase_inflates_rather_than_zeroes_the_orf():
    """ORF_CASE_SENSITIVITY: preserved from the C, unlike composition.

    No stop codon is recognised in lowercase, so the frame counter never
    resets and the whole contig reads as one open frame -- the maximum
    possible length, not zero.
    """
    rng = np.random.default_rng(7)
    seq = "".join(rng.choice(list("ACGT"), size=3000))

    upper = longest_orf(seq.upper())
    lower = longest_orf(seq.lower())

    assert lower > upper
    assert lower == 1000  # 3000 bases / 3, i.e. one unbroken frame


def test_lowercase_behaviour_matches_the_c_extension(oracle):
    rng = np.random.default_rng(11)
    seqs = []
    for _ in range(40):
        n = int(rng.integers(30, 900))
        s = "".join(rng.choice(list("ACGT"), size=n))
        seqs.extend([s.upper(), s.lower()])
    for seq, (_, want_orf) in zip(seqs, oracle(seqs)):
        assert longest_orf(seq) == want_orf, seq
