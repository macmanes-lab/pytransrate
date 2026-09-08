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
import mmap
import multiprocessing
import queue
from collections import OrderedDict

import numpy as np
import pysam

from pytransrate.assign import assign_fragments
from pytransrate.bam_metrics import (
    MalformedRecordStats,
    accumulate_into,
    accumulate_metrics,
    build_contigs,
    counts_buffer_size,
    estimate_realistic_distance,
    finalise_contigs,
    iter_alignments,
)
from pytransrate.segmenter import DEFAULT_NULL_PRIOR

__all__ = ["READ_STATS_KEYS", "ReadMetrics", "get_read_length"]

logger = logging.getLogger("pytransrate")

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

# ---------------------------------------------------------------------------
# STRIDING
#
# Assignment needs one fragment's records and the priors from quant.sf, and
# nothing else; fragments are contiguous in a read-ordered BAM. So fragments
# are independent units of work, and every accumulator they feed is additive
# -- coverage included, since it is carried as a difference array until it is
# integrated (COVERAGE_DIFF in bam_metrics). That makes the step separable.
#
# It is divided by striding: every worker reads the whole BAM and processes
# the fragments where `index % workers == worker`. That sounds wasteful, and
# it is -- each worker pays decompression and a fragment-boundary check on
# every record, ~0.26us of the ~1.56us a record costs to process. But the
# obvious alternative buys nothing. An unsorted BAM has no index, so exact
# byte ranges have to come from a scan that records bam.tell() at fragment
# boundaries, and that scan costs the same ~0.26us per record serially as
# striding's redundant decompression costs in parallel. Measured, the two
# come out level, and striding needs no seeking, no index and no guessing at
# record boundaries -- which matters in a file format where this module
# already documents an aligner writing records htslib refuses to parse.
#
# So the speedup is bounded by that fixed share, not by the worker count:
# ~3.7x at 8 workers, ~5.2x at 32, asymptotically ~6x. Beating it needs the
# parent to scan and dispatch ranges while workers run, which caps at ~6.1x
# for considerably more machinery, or splitting on BGZF block boundaries,
# which needs heuristic record-boundary detection. Neither is worth it.
#
# Workers write into anonymous mmaps created before the fork, which parent
# and children share; nothing large is pickled and nothing comes back through
# a pipe. That is also why this is fork-only, and why it falls back to the
# serial path rather than pretending elsewhere.
# ---------------------------------------------------------------------------

#: Per-contig accumulators the workers merge, in the order the shared scalar
#: buffer carries them.  Every one is a sum, which is what makes the merge a
#: single addition.  ``p_seq_true_sum`` is the only float; see PARALLEL_SUM.
_MERGED_FIELDS = (
    "reads_mapped",
    "bases_mapped",
    "fragments_mapped",
    "both_mapped",
    "properpair",
    "bridges",
    "good",
    "clipped_alignments",
    "clipped_bases",
    "leading_clipped_bases",
    "p_seq_true_sum",
)

#: The one accumulator above that is not an integer.
_FLOAT_FIELDS = frozenset({"p_seq_true_sum"})

# ---------------------------------------------------------------------------
# PARALLEL_SUM
#
# Every merged field is exact except p_seq_true_sum, which is a float sum
# over a contig's reads. Floating-point addition is not associative, so
# splitting a contig's reads across workers and adding the partials can
# differ in the last bits from adding them in BAM order.
#
# The merge is therefore done in worker order, which makes a run reproducible
# for a given --threads. It does not make it reproducible *across* thread
# counts: p_seq_true can move by ~1e-15 there, and it is a multiplicand of
# the contig score.
#
# That is fixable -- the per-read term simplifies to (35 - nm)/35, so the sum
# could be carried as two integers and divided once -- but doing so shifts
# today's numbers by ~1e-12, and this port has been deliberate about when it
# does that. Left as it is, and measured rather than assumed.
# ---------------------------------------------------------------------------

#: State handed to workers through the fork rather than through a pickle.
_WORKER: dict = {}

#: How long to wait on a worker's result before checking whether it is alive.
_RESULT_POLL_SECONDS = 5.0


def _worker_count(threads: int) -> int:
    """How many processes to actually use for ``--threads N``.

    Falls back to one wherever the fork-and-share scheme does not hold: the
    shared accumulators are anonymous mmaps inherited across ``fork``, which
    is Linux and macOS only.
    """
    if threads is None or threads < 2:
        return 1
    if "fork" not in multiprocessing.get_all_start_methods():
        logger.info("no fork available; assigning in one process")
        return 1
    return int(threads)


def _accumulate_stripe(worker_id: int) -> None:
    """Accumulate this worker's fragments into its shared buffers."""
    state = _WORKER
    references = state["references"]
    lengths = state["lengths"]
    bam_path = state["bam_path"]

    counts = np.frombuffer(state["counts"][worker_id], dtype=np.int32)
    scalars = np.frombuffer(state["scalars"][worker_id], dtype=np.float64)
    scalars = scalars.reshape(len(_MERGED_FIELDS), len(references))

    # Every worker reads every record, so they all meet the same malformed
    # ones; counting them once keeps the total honest. The run limit still
    # applies in each worker, so corruption still stops the run.
    stats = MalformedRecordStats() if worker_id == 0 else None
    error = None
    try:
        contigs = build_contigs(references, lengths, counts)
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            accumulate_into(
                contigs,
                assign_fragments(
                    iter_alignments(bam, bam_path, stats),
                    references,
                    state["expression"],
                    stride=state["workers"],
                    offset=worker_id,
                ),
                realistic_distance=state["realistic_distance"],
            )
        for row, field in enumerate(_MERGED_FIELDS):
            scalars[row] = [getattr(contig, field) for contig in contigs]
    except BaseException as exc:  # reported to the parent, then re-raised
        error = exc
    finally:
        state["results"].put(
            (worker_id, stats.skipped if stats else 0, _picklable(error))
        )
    if error is not None:
        raise error


def _picklable(error):
    """``error`` if it survives a pipe, else something that says what it was."""
    if error is None:
        return None
    try:
        import pickle

        pickle.loads(pickle.dumps(error))
    except Exception:
        return RuntimeError(f"{type(error).__name__}: {error}")
    return error


def _collect(processes, results_queue, bytes_per_worker: int) -> list:
    """One result per worker, without hanging on a worker that was killed.

    A worker that dies before it reports -- the OOM killer is the realistic
    way that happens, since the accumulators are sized by the assembly --
    would otherwise leave the parent blocked on a queue forever.  Waiting in
    slices and checking for a dead worker between them turns that into an
    error that says what to do about it.
    """
    results: list = []
    while len(results) < len(processes):
        try:
            results.append(results_queue.get(timeout=_RESULT_POLL_SECONDS))
        except queue.Empty:
            dead = [p for p in processes if p.exitcode not in (None, 0)]
            if dead:
                raise RuntimeError(
                    f"an assignment worker died with exit code {dead[0].exitcode} "
                    "before reporting. The usual cause is the machine running "
                    f"out of memory, at {bytes_per_worker / 1e9:.1f} GB per "
                    "worker -- lower --threads"
                ) from None
    return results


def _accumulate_parallel(
    bam_path,
    references,
    lengths,
    expression,
    realistic_distance: int,
    nullprior: float,
    workers: int,
    malformed: MalformedRecordStats,
):
    """Assign and accumulate across ``workers`` processes. See STRIDING."""
    context = multiprocessing.get_context("fork")
    n_contigs = len(references)
    counts_bytes = counts_buffer_size(lengths) * 4
    scalar_bytes = len(_MERGED_FIELDS) * n_contigs * 8

    logger.info(
        "assigning across %d processes (%.1f GB of shared accumulators)",
        workers,
        workers * (counts_bytes + scalar_bytes) / 1e9,
    )

    counts_blocks = [mmap.mmap(-1, counts_bytes) for _ in range(workers)]
    scalar_blocks = [mmap.mmap(-1, scalar_bytes) for _ in range(workers)]
    results_queue = context.Queue()
    processes = []
    try:
        _WORKER.update(
            bam_path=str(bam_path),
            references=references,
            lengths=lengths,
            expression=expression,
            realistic_distance=realistic_distance,
            workers=workers,
            counts=counts_blocks,
            scalars=scalar_blocks,
            results=results_queue,
        )
        processes = [
            context.Process(target=_accumulate_stripe, args=(worker_id,))
            for worker_id in range(workers)
        ]
        for process in processes:
            process.start()
        results = _collect(processes, results_queue, counts_bytes + scalar_bytes)
        for process in processes:
            process.join()
    finally:
        _WORKER.clear()
        for process in processes:
            if process.is_alive():
                process.terminate()
        results_queue.close()

    for _worker_id, skipped, error in results:
        if error is not None:
            raise error
        malformed.skipped += skipped

    # Sum into the first worker's buffers, then take a private copy: the
    # contigs below hold views into it, and the mmaps do not outlive this
    # function. Peak cost is one extra buffer.
    merged_counts = np.frombuffer(counts_blocks[0], dtype=np.int32)
    for block in counts_blocks[1:]:
        merged_counts += np.frombuffer(block, dtype=np.int32)
    merged_counts = merged_counts.copy()

    shape = (len(_MERGED_FIELDS), n_contigs)
    merged_scalars = np.frombuffer(scalar_blocks[0], dtype=np.float64).reshape(shape)
    for block in scalar_blocks[1:]:
        merged_scalars += np.frombuffer(block, dtype=np.float64).reshape(shape)
    merged_scalars = merged_scalars.copy()

    for block in counts_blocks + scalar_blocks:
        block.close()

    contigs = build_contigs(references, lengths, merged_counts)
    for row, field in enumerate(_MERGED_FIELDS):
        values = merged_scalars[row]
        if field not in _FLOAT_FIELDS:
            values = values.astype(np.int64)
        for contig, value in zip(contigs, values.tolist()):
            setattr(contig, field, value)

    finalise_contigs(contigs, nullprior=nullprior)
    return contigs

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
        threads: int = 1,
    ) -> "ReadMetrics":
        """Assign fragments, accumulate per-contig metrics, and aggregate.

        Args:
            bam_path: read-ordered BAM from the aligner, still carrying every
                multi-mapping alignment.
            expression: parsed ``quant.sf``.
            fragments: total fragments in the library, from the aligner.
            read_length: max read length, used for the coverage estimate.
            nullprior: prior handed to the segmenter.
            threads: processes to divide the fragments across.  See STRIDING
                for what this does and does not buy; 1 keeps everything in
                this process.
        """
        self.fragments = fragments
        if read_length:
            self.read_length = read_length

        realistic_distance = estimate_realistic_distance(str(bam_path))
        logger.debug("realistic fragment distance: %d", realistic_distance)

        malformed = MalformedRecordStats()
        workers = _worker_count(threads)
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            references = list(bam.references)
            lengths = list(bam.lengths)

        if workers > 1:
            contig_metrics = _accumulate_parallel(
                bam_path,
                references,
                lengths,
                expression,
                realistic_distance=realistic_distance,
                nullprior=nullprior,
                workers=workers,
                malformed=malformed,
            )
        else:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                assigned = assign_fragments(
                    iter_alignments(bam, str(bam_path), malformed),
                    references,
                    expression,
                )
                contig_metrics = accumulate_metrics(
                    references,
                    lengths,
                    assigned,
                    realistic_distance=realistic_distance,
                    nullprior=nullprior,
                )

        if malformed.skipped:
            logger.warning("%s", malformed.describe())

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
