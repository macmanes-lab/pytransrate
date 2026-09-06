"""Assembly scoring and cutoff optimisation.

Port of ``lib/transrate/score_optimiser.rb``.

The optimiser walks contigs shortest-score-first, removing each in turn, and
records the assembly score that would result from cutting there.  The cutoff
maximising that score is the "optimal" one, and is what
``classify_contigs`` then uses to split good from bad.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass

__all__ = ["ScoreOptimiser", "geomean"]


def geomean(values) -> float:
    """Geometric mean, via logs. ``Contig#score`` is floored at 0.01."""
    values = list(values)
    if not values:
        return 0.0
    total = 0.0
    for value in values:
        total += math.log(value)
    return math.exp(total / len(values))


@dataclass
class ScoreOptimiser:
    """Computes the raw, weighted and cutoff-optimised assembly scores.

    Args:
        assembly: an :class:`~transrate.assembly.Assembly`.
        fragments: total read fragments (``read_stats['fragments']``).
        good: fragments mapping consistently (``read_stats['good_mappings']``).
    """

    assembly: object
    fragments: int
    good: int

    _optimal: float | None = None
    _cutoff: float | None = None

    def raw_score(self) -> float:
        """Geometric mean contig score, scaled by the good-mapping rate."""
        contigs = list(self.assembly.values())
        if not contigs or not self.fragments:
            return 0.0
        contig_score = geomean(c.score for c in contigs)
        return contig_score * (self.good / self.fragments)

    def weighted_score(self) -> float:
        """As :meth:`raw_score` but weighting each contig by expression."""
        contigs = list(self.assembly.values())
        if not contigs or not self.fragments:
            return 0.0
        total = sum(c.score * c.tpm for c in contigs)
        contig_weighted = total / len(contigs)
        return contig_weighted * (self.good / self.fragments)

    def optimal_score(self, csv_path: str | None = None) -> tuple[float, float]:
        """Find the cutoff maximising the assembly score.

        .. note:: **CUTOFF_BOUNDARY.** The returned cutoff is the score of the
           last contig *removed* to reach the optimum, but
           :meth:`Contig.classify` keeps contigs with ``score >= cutoff``.
           Contigs sitting exactly on the cutoff are therefore classified
           good even though the optimum was computed with them excluded.
           This bites hardest at the 0.01 score floor, where many dead
           contigs share a score and collapse to a single cutoff key: they
           all survive classification.  Reproduced from the Ruby, which
           behaves identically -- ORP's good/bad split has always had this
           boundary.

        Args:
            csv_path: where to write the cutoff/score curve.  The Ruby always
                wrote ``<prefix>_score_optimisation.csv``; omit to skip.

        Returns:
            ``(optimal_score, cutoff)``.
        """
        if self._optimal is not None:
            return self._optimal, self._cutoff

        contigs = list(self.assembly.values())
        if not contigs or not self.fragments:
            self._optimal, self._cutoff = 0.0, 0.0
            return self._optimal, self._cutoff

        product = 0.0
        good = 0
        for contig in contigs:
            product += math.log(contig.score)
            good += contig.good
        count = len(contigs)

        # Keyed by contig score: contigs sharing a score collapse to one
        # cutoff, last one winning, as the Ruby Hash did.
        cutoff_scores: dict[float, float] = {}
        for contig in sorted(contigs, key=lambda c: c.score):
            product -= math.log(contig.score)
            good -= contig.good
            count -= 1
            if count <= 0:
                # The Ruby divided by zero here and got Infinity/NaN, which
                # then lost every `score > optimal` comparison. Skipping is
                # equivalent and keeps the arithmetic defined.
                continue
            score = math.exp(product / count) * (good / self.fragments)
            cutoff_scores[contig.score] = score

        optimal = 0.0
        cutoff = 0.0
        rows = []
        for contig_score, score in cutoff_scores.items():
            rows.append((round(contig_score, 5), round(score, 5)))
            if score > optimal:
                optimal = score
                cutoff = contig_score

        if csv_path is not None:
            with open(csv_path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["cutoff", "assembly_score"])
                writer.writerows(rows)

        self._optimal, self._cutoff = optimal, cutoff
        return optimal, cutoff
