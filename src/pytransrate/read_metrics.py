"""Aggregate read-based metrics for an assembly.

Port of ``lib/transrate/read_metrics.rb``.  The Ruby drove three external
steps here -- snap, salmon with ``--sampleOut``, then the ``bam-read``
binary over the resulting ``postSample.bam``.  Only the first two are still
external; assignment and per-contig accumulation happen in-process via
:mod:`~pytransrate.assign` and :mod:`~pytransrate.bam_metrics`.
"""

from __future__ import annotations

import gzip
import logging
from collections import OrderedDict

import pysam

from pytransrate.assign import assign_fragments
from pytransrate.bam_metrics import (
    accumulate_metrics,
    estimate_realistic_distance,
    iter_alignments,
)
from pytransrate.segmenter import DEFAULT_NULL_PRIOR

__all__ = ["READ_STATS_KEYS", "ReadMetrics", "get_read_length"]

logger = logging.getLogger("transrate")

#: Key order of ``ReadMetrics#read_stats``; mirrored in pytransrate.output.
READ_STATS_KEYS = (
    "fragments",
    "fragments_mapped",
    "p_fragments_mapped",
    "good_mappings",
    "p_good_mapping",
    "bad_mappings",
    "potential_bridges",
    "bases_uncovered",
    "p_bases_uncovered",
    "contigs_uncovbase",
    "p_contigs_uncovbase",
    "contigs_uncovered",
    "p_contigs_uncovered",
    "contigs_lowcovered",
    "p_contigs_lowcovered",
    "contigs_segmented",
    "p_contigs_segmented",
)

#: Reads inspected when estimating the maximum read length.
_READ_LENGTH_SAMPLE = 5000

#: Mean coverage below which a contig counts as uncovered / low-covered.
_UNCOVERED_BELOW = 1
_LOWCOVERED_BELOW = 10

#: A contig is segmented when p_not_segmented drops below this.
_SEGMENTED_BELOW = 0.5


def _open_maybe_gzip(path):
    if str(path).endswith((".gz", ".gzip")):
        return gzip.open(path, "rt")
    return open(path)


def get_read_length(reads) -> int:
    """Maximum read length over the first few thousand reads.

    Matches ``ReadMetrics#get_read_length``, which sampled the first file
    only.  Unlike the Ruby this also reads gzipped FASTQ.
    """
    first = str(reads).split(",")[0]
    longest = 0
    with _open_maybe_gzip(first) as handle:
        for index, line in enumerate(handle):
            if index >= _READ_LENGTH_SAMPLE * 4:
                break
            if index % 4 == 1:
                longest = max(longest, len(line.strip()))
    return longest


class ReadMetrics:
    """Read-mapping metrics for one assembly."""

    def __init__(self, assembly):
        self.assembly = assembly
        self.has_run = False
        self.read_length = 100

        self.fragments = 0
        self.fragments_mapped = 0
        self.good = 0
        self.bad = 0
        self.potential_bridges = 0
        self.bases_uncovered = 0
        self.contigs_uncovbase = 0
        self.contigs_uncovered = 0
        self.contigs_lowcovered = 0
        self.contigs_segmented = 0

        self.alignments = 0
        self.clipped_alignments = 0
        self.clipped_bases = 0
        self.leading_clipped_bases = 0

    # -- driving ----------------------------------------------------------

    def run(
        self,
        bam_path,
        expression,
        fragments: int,
        read_length: int | None = None,
        nullprior: float = DEFAULT_NULL_PRIOR,
    ) -> "ReadMetrics":
        """Assign fragments, accumulate per-contig metrics, and aggregate.

        Args:
            bam_path: read-ordered BAM from the aligner, still carrying every
                multi-mapping alignment.
            expression: parsed ``quant.sf``.
            fragments: total fragments in the library, from the aligner.
            read_length: max read length, used for the coverage estimate.
            nullprior: prior handed to the segmenter.
        """
        self.fragments = fragments
        if read_length:
            self.read_length = read_length

        realistic_distance = estimate_realistic_distance(str(bam_path))
        logger.debug("realistic fragment distance: %d", realistic_distance)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            references = list(bam.references)
            assigned = assign_fragments(
                iter_alignments(bam, str(bam_path)), references, expression
            )
            contig_metrics = accumulate_metrics(
                references,
                bam.lengths,
                assigned,
                realistic_distance=realistic_distance,
                nullprior=nullprior,
            )

        self._summarise_clipping(contig_metrics)
        self._populate_contigs(contig_metrics)
        self._analyse_expression(expression)
        self._update_proportions()
        self.has_run = True
        return self

    # -- aggregation ------------------------------------------------------

    def _summarise_clipping(self, contig_metrics) -> None:
        """Report soft clipping, which the BAM would otherwise take with it.

        snap 2.x clips reads that hang over a contig end, so those terminal
        bases lose their support and sCcov falls. That is real, but it is
        invisible once the BAM is deleted, and it is the first thing anyone
        comparing scores against the Ruby needs to know.
        """
        for metrics in contig_metrics:
            self.alignments += metrics.reads_mapped
            self.clipped_alignments += metrics.clipped_alignments
            self.clipped_bases += metrics.clipped_bases
            self.leading_clipped_bases += metrics.leading_clipped_bases

        if not self.alignments:
            return

        fraction = self.clipped_alignments / self.alignments
        logger.info(
            "soft-clipped alignments: %d / %d (%.1f%%), %d bases clipped",
            self.clipped_alignments,
            self.alignments,
            100 * fraction,
            self.clipped_bases,
        )
        if self.clipped_alignments:
            logger.info(
                "mean soft clip: %.1f bp (%.1f bp leading)",
                self.clipped_bases / self.clipped_alignments,
                self.leading_clipped_bases / self.clipped_alignments,
            )

    def _populate_contigs(self, contig_metrics) -> None:
        """Fold per-contig BAM metrics into the assembly's Contig objects."""
        for metrics in contig_metrics:
            if metrics.name not in self.assembly:
                continue
            contig = self.assembly[metrics.name]

            contig.p_seq_true = metrics.p_seq_true()
            contig.set_uncovered_bases(metrics.bases_uncovered)
            self.bases_uncovered += metrics.bases_uncovered

            # The Ruby guarded on >1, so single-fragment contigs keep p_good 0.
            if metrics.fragments_mapped and metrics.fragments_mapped > 1:
                contig.p_good = metrics.good / metrics.fragments_mapped

            contig.p_not_segmented = metrics.p_not_segmented
            if contig.p_not_segmented < _SEGMENTED_BELOW:
                self.contigs_segmented += 1

            contig.in_bridges = metrics.bridges
            if metrics.bridges > 1:
                self.potential_bridges += 1

            self.fragments_mapped += metrics.fragments_mapped
            contig.good = metrics.good
            self.good += metrics.good

            if metrics.bases_uncovered > 0:
                self.contigs_uncovbase += 1

        self.bad = self.fragments_mapped - self.good

    def _analyse_expression(self, expression) -> None:
        """Attach salmon's estimates and derive mean coverage."""
        for name, values in (expression or {}).items():
            contig_name = name.split()[0].split("|")[0].rstrip(";")
            if contig_name not in self.assembly:
                continue
            contig = self.assembly[contig_name]

            eff_len = values["eff_len"]
            if eff_len == 0:
                coverage = 0.0
            else:
                coverage = values["eff_count"] * self.read_length / eff_len

            if coverage < _UNCOVERED_BELOW:
                self.contigs_uncovered += 1
            if coverage < _LOWCOVERED_BELOW:
                self.contigs_lowcovered += 1

            contig.coverage = round(coverage, 2)
            contig.eff_length = eff_len
            contig.eff_count = values["eff_count"]
            contig.tpm = values["tpm"]

    def _update_proportions(self) -> None:
        n_bases = float(self.assembly.n_bases) or 1.0
        n_contigs = float(self.assembly.size) or 1.0
        fragments = float(self.fragments) or 1.0

        self.p_bases_uncovered = self.bases_uncovered / n_bases
        self.p_contigs_uncovbase = self.contigs_uncovbase / n_contigs
        self.p_contigs_uncovered = self.contigs_uncovered / n_contigs
        self.p_contigs_lowcovered = self.contigs_lowcovered / n_contigs
        self.p_contigs_segmented = self.contigs_segmented / n_contigs
        self.p_good_mapping = self.good / fragments
        self.p_fragments_mapped = self.fragments_mapped / fragments

    # -- output -----------------------------------------------------------

    def read_stats(self) -> "OrderedDict[str, float]":
        stats: "OrderedDict[str, float]" = OrderedDict()
        stats["fragments"] = self.fragments
        stats["fragments_mapped"] = self.fragments_mapped
        stats["p_fragments_mapped"] = self.p_fragments_mapped
        stats["good_mappings"] = self.good
        stats["p_good_mapping"] = self.p_good_mapping
        stats["bad_mappings"] = self.bad
        stats["potential_bridges"] = self.potential_bridges
        stats["bases_uncovered"] = self.bases_uncovered
        stats["p_bases_uncovered"] = self.p_bases_uncovered
        stats["contigs_uncovbase"] = self.contigs_uncovbase
        stats["p_contigs_uncovbase"] = self.p_contigs_uncovbase
        stats["contigs_uncovered"] = self.contigs_uncovered
        stats["p_contigs_uncovered"] = self.p_contigs_uncovered
        stats["contigs_lowcovered"] = self.contigs_lowcovered
        stats["p_contigs_lowcovered"] = self.p_contigs_lowcovered
        stats["contigs_segmented"] = self.contigs_segmented
        stats["p_contigs_segmented"] = self.p_contigs_segmented
        return stats
