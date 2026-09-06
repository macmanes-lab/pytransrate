"""A single contig and its metrics.

Port of ``lib/transrate/contig.rb``.  The score components and their 0.01
floor are reproduced exactly -- they define the transrate score, and ORP
selects contigs on it via ``scripts/pick_best_contigs.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pytransrate.sequence import base_composition, longest_orf

__all__ = ["SCORE_FLOOR", "Contig"]

#: Every score component is floored here before multiplying, so a single
#: zero component cannot collapse the product to zero.
SCORE_FLOOR = 0.01


@dataclass
class Contig:
    """One assembled sequence plus whatever metrics have been computed."""

    name: str
    seq: str

    # Read-based metrics, filled in from the BAM and salmon.
    eff_length: float = 0.0
    eff_count: float = 0.0
    tpm: float = 0.0
    coverage: float = 0.0
    uncovered_bases: int = 0
    p_uncovered_bases: float = 1.0
    p_seq_true: float = 0.0
    in_bridges: int = 0
    p_good: float = 0.0
    p_not_segmented: float = 1.0
    good: int = 0
    classification: str = "unknown"

    # Reference-based metrics (unused unless --reference is given).
    has_crb: bool = False
    reference_coverage: float = 0.0
    hits: list = field(default_factory=list)

    _base_composition: dict | None = field(default=None, repr=False)
    _orf_length: int | None = field(default=None, repr=False)

    def __post_init__(self):
        # The Ruby strips NULs from the sequence and trailing semicolons from
        # the identifier, because BLAST strips the latter.
        self.seq = self.seq.replace("\0", "")
        self.name = self.name.rstrip(";")
        if self.uncovered_bases == 0 and len(self.seq):
            # Matches Contig#initialize: assume fully uncovered until told
            # otherwise.
            self.uncovered_bases = len(self.seq)
            self.p_uncovered_bases = 1.0

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def length(self) -> int:
        return len(self.seq)

    # -- composition ------------------------------------------------------

    @property
    def composition(self) -> dict[str, int]:
        if self._base_composition is None:
            self._base_composition = base_composition(self.seq)
        return self._base_composition

    @property
    def bases_a(self) -> int:
        return self.composition["a"]

    @property
    def bases_c(self) -> int:
        return self.composition["c"]

    @property
    def bases_g(self) -> int:
        return self.composition["g"]

    @property
    def bases_t(self) -> int:
        return self.composition["t"]

    @property
    def bases_n(self) -> int:
        return self.composition["n"]

    @property
    def prop_gc(self) -> float:
        if not self.length:
            return 0.0
        return (self.bases_g + self.bases_c) / self.length

    @property
    def orf_length(self) -> int:
        """Longest ORF, in codons. See ORF_CASE_SENSITIVITY."""
        if self._orf_length is None:
            self._orf_length = longest_orf(self.seq)
        return self._orf_length

    # -- coverage ---------------------------------------------------------

    def set_uncovered_bases(self, n: int) -> None:
        self.uncovered_bases = n
        self.p_uncovered_bases = n / self.length if self.length else 1.0

    @property
    def p_bases_covered(self) -> float:
        return 1.0 - self.p_uncovered_bases

    # -- score ------------------------------------------------------------

    @property
    def score(self) -> float:
        """Product of the four score components, each floored at 0.01.

        These are sCcov, sCseg, sCord and sCnuc in the transrate paper.
        """
        product = (
            max(self.p_bases_covered, SCORE_FLOOR)
            * max(self.p_not_segmented, SCORE_FLOOR)
            * max(self.p_good, SCORE_FLOOR)
            * max(self.p_seq_true, SCORE_FLOOR)
        )
        return max(product, SCORE_FLOOR)

    def alt_scores(self) -> dict[str, float]:
        """Leave-one-out scores, one per component.

        ``Contig#alt_score`` in the Ruby; used by ``classify_contigs`` to
        write the ``single_component_bad`` breakdown.
        """
        parts = {
            "cov": max(self.p_bases_covered, SCORE_FLOOR),
            "seg": max(self.p_not_segmented, SCORE_FLOOR),
            "good": max(self.p_good, SCORE_FLOOR),
            "seq": max(self.p_seq_true, SCORE_FLOOR),
        }
        out = {}
        for missing in parts:
            product = 1.0
            for key, value in parts.items():
                if key != missing:
                    product *= value
            out[missing] = product
        return out

    def classify(self, cutoff: float) -> str:
        self.classification = "good" if self.score >= cutoff else "bad"
        return self.classification

    # -- output rows ------------------------------------------------------

    def basic_metrics(self) -> dict:
        return {
            "length": self.length,
            "prop_gc": self.prop_gc,
            "orf_length": self.orf_length,
        }

    def read_metrics(self) -> dict:
        return {
            "in_bridges": self.in_bridges,
            "p_good": self.p_good,
            "p_bases_covered": self.p_bases_covered,
            "p_seq_true": self.p_seq_true,
            "score": self.score,
            "p_not_segmented": self.p_not_segmented,
            "eff_length": self.eff_length,
            "eff_count": self.eff_count,
            "tpm": self.tpm,
            "coverage": self.coverage,
            "sCnuc": self.p_seq_true,
            "sCcov": self.p_bases_covered,
            "sCord": self.p_good,
            "sCseg": self.p_not_segmented,
        }

    def comparative_metrics(self) -> dict:
        if self.has_crb:
            return {
                "has_crb": True,
                "reference_coverage": self.reference_coverage,
                "hits": ";".join(str(h) for h in self.hits),
            }
        return {"has_crb": False, "reference_coverage": "NA", "hits": "NA"}

    def to_fasta(self) -> str:
        return f">{self.name}\n{self.seq}\n"
