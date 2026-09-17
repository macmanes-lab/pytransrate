"""Command-line interface.

Port of ``lib/transrate/cmdline.rb``.  The option names ORP depends on are
unchanged -- ``-a/--assembly``, ``-o/--output``, ``-t/--threads``,
``--left``, ``--right`` -- so existing invocations keep working:

    pytransrate -o <dir> -t <cpu> -a <fasta> --left <r1> --right <r2>

``--install-deps`` is gone.  It drove the ``bindeps`` gem, which fetched
binaries from URLs that no longer resolve; snap-aligner and salmon are now
ordinary conda/bioconda packages.

``--reference`` is not implemented in this port.  The Ruby routed it through
the unmaintained ``crb-blast`` gem, and ORP never passes it.  The flag is
still parsed so that a script using it gets a clear error rather than
silently different output.
"""

from __future__ import annotations

import argparse
import logging
import os
import shlex
import sys
from pathlib import Path

from pytransrate import __version__
from pytransrate.assembly import Assembly, AssemblyError
from pytransrate.bam_metrics import MalformedBamError
from pytransrate.banner import TAGLINE, print_banner
from pytransrate.cmd import CommandError
from pytransrate.compression import plain_path, strip_gzip_suffix
from pytransrate.mapper import Snap, align_done_path
from pytransrate.memory import parse_size
from pytransrate.output import (
    READ_STATS_KEYS,
    write_assemblies_csv,
    write_contigs_csv,
)
from pytransrate.quantify import Salmon
from pytransrate.read_metrics import ReadMetrics, get_read_length
from pytransrate.score import ScoreOptimiser

logger = logging.getLogger("pytransrate")

#: Written into the output directory, as the Ruby did.
ASSEMBLIES_CSV = "assemblies.csv"
CONTIGS_CSV = "contigs.csv"

# ---------------------------------------------------------------------------
# REPORTED_METRICS
#
# Everything computed goes to assemblies.csv; these are the ones worth
# reading while a run is in progress or out of a log afterwards. The contig
# block is deliberately three lines -- min, max and N50 -- because the rest
# of the length distribution is not what anyone watches a run for, and the
# CSV has it. The mapping block is complete: it is what changes when the
# aligner settings change, so leaving any of it out would mean going to the
# CSV to answer the obvious follow-up question.
#
# Values are formatted exactly as write_assemblies_csv rounds them (5 places),
# so a number read off the log matches the one in the CSV.
# ---------------------------------------------------------------------------

#: Contig statistics reported at INFO. Keys of ``Assembly.basic_stats()``.
REPORTED_CONTIG_KEYS = (
    ("smallest", "min contig length"),
    ("largest", "max contig length"),
    ("n50", "N50"),
)


_EXAMPLES = """
examples:
  # sequence metrics only -- no reads, no aligner needed
  pytransrate -a assembly.fa -o results

  # the full analysis: contig, read-mapping and score metrics
  pytransrate -a assembly.fa --left r1.fq --right r2.fq -t 16 -o results

  # several assemblies, each into its own subdirectory of results/
  pytransrate -a one.fa,two.fa --left r1.fq --right r2.fq -o results

  # multiple read files, given as matched comma-separated lists
  pytransrate -a assembly.fa --left a1.fq,b1.fq --right a2.fq,b2.fq -o results

  # show every external command before it runs
  pytransrate -a assembly.fa --left r1.fq --right r2.fq -o results --loglevel debug

output:
  assemblies.csv   assembly-level metrics; score and optimal_score are the
                   37th and 38th columns
  contigs.csv      per-contig metrics; score is the 9th column
  *_score_optimisation.csv
                   the cutoff/score curve the optimiser walked

Column order in both files is an interface, not presentation -- downstream
tools read them positionally.
"""


class _HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Wider option column, so the grouped flags stay readable."""

    def __init__(self, prog):
        super().__init__(prog, max_help_position=32)


def _memory_size(text: str) -> int:
    """``--max-memory`` as bytes, with argparse's own error on nonsense."""
    try:
        return parse_size(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pytransrate",
        formatter_class=_HelpFormatter,
        description=(
            f"pytransrate: {TAGLINE}.\n\n"
            "Maps the reads back to the assembly, quantifies expression, and\n"
            "reduces the result to a score per contig and one for the assembly,\n"
            "measuring how well the assembly is supported by its own reads."
        ),
        epilog=_EXAMPLES,
    )

    required = parser.add_argument_group("required")
    required.add_argument(
        "-a", "--assembly",
        required=True,
        metavar="FASTA",
        help=(
            "assembly file(s) in FASTA format, comma-separated. "
            "gzipped files are read directly"
        ),
    )

    reads = parser.add_argument_group(
        "reads",
        "Give both to enable read-mapping metrics and the transrate score.\n"
        "Without them, only sequence-based metrics are computed.",
    )
    reads.add_argument(
        "--left",
        metavar="FASTQ",
        help="left reads, comma-separated; gzipped files are read directly",
    )
    reads.add_argument(
        "--right", metavar="FASTQ", help="right reads, comma-separated"
    )

    general = parser.add_argument_group("general")
    general.add_argument(
        "-o", "--output",
        default="transrate_results",
        metavar="DIR",
        help="output directory (default: transrate_results)",
    )
    general.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        metavar="N",
        help="threads to use (default: 8)",
    )
    general.add_argument(
        # --mem because that is what a scheduler directive and the pipelines
        # wrapping this call it; a job that already knows it asked for 670 GB
        # should be able to pass that figure straight through.
        "--max-memory", "--mem",
        type=_memory_size,
        default=None,
        metavar="SIZE",
        help="memory the read-metrics step may use, e.g. 670G (a bare number "
             "is GB, so --mem 670 is 670 GB). Caps the processes it forks, "
             "since each one holds a full-size copy of the per-base coverage "
             "accumulators. The default reads the cgroup, the Slurm "
             "allocation and /proc/meminfo, which is right unless your "
             "scheduler or pipeline enforces a limit none of those show",
    )
    general.add_argument(
        "--loglevel",
        default="info",
        choices=["error", "warn", "info", "debug"],
        help="logging verbosity (default: info); debug logs every command",
    )
    general.add_argument(
        "--keep-bam",
        action="store_true",
        help="keep the alignment BAM instead of deleting it on success",
    )
    general.add_argument(
        "--no-banner", action="store_true", help="suppress the startup banner"
    )
    general.add_argument(
        "--version", action="version", version=f"pytransrate {__version__}"
    )

    index = parser.add_argument_group(
        "snap index tuning",
        "Only needed when a large or repetitive assembly exhausts\nsnap's genome locations.",
    )
    index.add_argument(
        "--location-size",
        type=int,
        choices=range(4, 9),
        metavar="{4-8}",
        default=None,
        help=(
            "snap -locationSize. Default sweeps 4 up to 8 when snap runs "
            "out of genome locations; each failed attempt is a full index "
            "build, so set this if you already know the value"
        ),
    )
    index.add_argument(
        "--seed-size",
        type=int,
        default=23,
        metavar="N",
        help=(
            "snap index -s (default: 23). The other fix when an assembly "
            "is still too big at every location size"
        ),
    )
    index.add_argument(
        "--padding",
        type=int,
        default=1000,
        metavar="N",
        help=(
            "snap index -p, Ns between contigs (default: 1000; snap's own "
            "default is 2000). Counted toward the genome size, so on an "
            "assembly with millions of contigs the padding, not the "
            "sequence, is what exhausts the location namespace. Must stay "
            "above --edit-distance"
        ),
    )

    mapping = parser.add_argument_group(
        "snap mapping tuning",
        "Defaults are the configuration verified against real assemblies.\n"
        "The original Ruby's values crash snap 2.x with SIGFPE; see\n"
        "MULTI_ALIGNMENT_SETTINGS in mapper.py.",
    )
    mapping.add_argument(
        "--multi-edit-distance",
        type=int,
        default=2,
        metavar="N",
        help=(
            "snap -om (default: 2): extra edit distance admitted for "
            "secondary alignments, which is what fragment assignment chooses "
            "between. Higher values explode on a redundant assembly for "
            "almost no extra signal"
        ),
    )
    mapping.add_argument(
        "--max-alignments-per-contig",
        type=int,
        default=1,
        metavar="N",
        help=(
            "snap -mpc (default: 1), applied before -omax: the best placement "
            "per candidate contig. 0 disables the cap. Raising it inflates "
            "the score without improving the assembly -- measured on three "
            "assemblies, see MULTI_ALIGNMENT_SETTINGS in mapper.py"
        ),
    )
    mapping.add_argument(
        "--max-alignments-per-pair",
        type=int,
        default=10,
        metavar="N",
        help="snap -omax, cap on alignments per pair (default: 10)",
    )
    mapping.add_argument(
        "--extra-search-depth",
        type=int,
        default=2,
        metavar="N",
        help="snap -D (default: 2). Must be >= --multi-edit-distance",
    )
    mapping.add_argument(
        "--max-seed-hits",
        type=int,
        default=4000,
        metavar="N",
        help="snap -H (default: 4000, snap's own default; the Ruby used 300000)",
    )
    mapping.add_argument(
        "--edit-distance",
        type=int,
        default=30,
        metavar="N",
        help="snap -d, max edit distance per pair (default: 30)",
    )
    mapping.add_argument(
        "--max-candidate-pool",
        type=int,
        default=None,
        metavar="N",
        help=(
            "snap -mcp. Not passed by default: the Ruby's value overflowed "
            "snap's atoi() into an arbitrary number. Must be under 2147483647"
        ),
    )

    quant = parser.add_argument_group("salmon")
    quant.add_argument(
        "--no-error-model",
        action="store_true",
        help=(
            "don't pass --errorModel to salmon. Only sensible if your BAM "
            "carries AS tags; snap-aligner does not emit them"
        ),
    )

    unsupported = parser.add_argument_group("not implemented")
    unsupported.add_argument(
        "-r", "--reference",
        metavar="FASTA",
        help=(
            "reference proteome/transcriptome. Not implemented in this port; "
            "passing it raises rather than silently changing the output"
        ),
    )
    return parser


def configure_logging(level: str) -> None:
    """Send the whole run report to stdout.

    Diagnostics and results share one stream on purpose: they are read
    together, out of a single redirected log, and splitting them across
    stdout and stderr interleaves them unpredictably in that file.
    """
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)5s] %(asctime)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )


def format_metric(value) -> str:
    """Format one metric as write_assemblies_csv would round it."""
    if isinstance(value, bool) or value is None:
        return str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.5f}"
    return str(value)


def log_metrics(title: str, values: dict, keys) -> None:
    """Report a labelled block of metrics at INFO. See REPORTED_METRICS."""
    pairs = []
    for key in keys:
        name, label = key if isinstance(key, tuple) else (key, key)
        if name in values:
            pairs.append((label, format_metric(values[name])))
    if not pairs:
        return

    width = max(len(label) for label, _ in pairs)
    logger.info("%s:", title)
    for label, value in pairs:
        logger.info("  %-*s  %s", width, label, value)


def check_arguments(args) -> list[str]:
    """Validate inputs, returning the expanded assembly paths."""
    if args.reference:
        raise CommandError(
            "--reference is not implemented in this port. The Ruby "
            "implementation routed it through the unmaintained crb-blast "
            "gem; reference-based metrics are not available here."
        )

    assemblies = []
    for path in args.assembly.split(","):
        full = os.path.abspath(os.path.expanduser(path))
        if not os.path.exists(full):
            raise CommandError(f"assembly fasta file does not exist: {full}")
        assemblies.append(full)

    if len(set(assemblies)) != len(assemblies):
        raise CommandError(
            "the same assembly was supplied more than once to --assembly; "
            "output paths are derived from it, so each must appear once"
        )

    if bool(args.left) != bool(args.right):
        raise CommandError("--left and --right must be given together")

    if args.left:
        left = args.left.split(",")
        right = args.right.split(",")
        if len(left) != len(right):
            raise CommandError(
                "please provide the same number of left and right read files"
            )
        for path in left + right:
            if not os.path.exists(os.path.expanduser(path)):
                raise CommandError(f"read fastq file does not exist: {path}")
        # Resolve now: analysis runs with the cwd changed into the result
        # directory, so relative paths would resolve against the wrong root.
        args.left = ",".join(
            os.path.abspath(os.path.expanduser(p)) for p in left
        )
        args.right = ",".join(
            os.path.abspath(os.path.expanduser(p)) for p in right
        )

    return assemblies


def analyse_assembly(
    assembly_path, args, result_dir: Path, defer_delete: list | None = None
) -> dict:
    """Run every requested metric for one assembly."""
    logger.info("loading assembly: %s", assembly_path)
    assembly = Assembly(assembly_path)

    result = {"assembly": str(assembly_path)}
    logger.info("calculating contig metrics...")
    result.update(assembly.basic_stats())
    result.update(assembly.contig_metrics())
    log_metrics("contig metrics", result, REPORTED_CONTIG_KEYS)

    with_reads = bool(args.left and args.right)
    if not with_reads:
        logger.info("no reads provided, skipping read diagnostics")
        write_contigs_csv(
            assembly, str(result_dir / CONTIGS_CSV), with_reads=False
        )
        return result

    left, right = args.left, args.right  # already absolute, see check_arguments

    # snap-aligner and salmon are given a filename, not a handle, so a
    # gzipped assembly is decompressed for the length of this block and
    # dropped again. Assembly above read the gzip itself, and without reads
    # neither tool runs, so nothing is decompressed for a sequence-only run.
    # The reads stay gzipped: snap reads them that way, and a library is
    # orders of magnitude larger than an assembly. See pytransrate.compression.
    snap = Snap()
    with plain_path(assembly_path, result_dir) as fasta:
        logger.info("mapping reads with snap-aligner...")
        snap.build_index(
            fasta,
            threads=args.threads,
            seed_size=args.seed_size,
            location_size=args.location_size,
            padding=args.padding,
        )
        bam = snap.map_reads(
            left, right, threads=args.threads,
            max_seed_hits=args.max_seed_hits,
            edit_distance=args.edit_distance,
            extra_search_depth=args.extra_search_depth,
            multi_edit_distance=args.multi_edit_distance,
            max_alignments_per_pair=args.max_alignments_per_pair,
            max_alignments_per_contig=(
                args.max_alignments_per_contig
                if args.max_alignments_per_contig > 0
                else None
            ),
            max_candidate_pool=args.max_candidate_pool,
        )
        logger.info("%d fragments in library", snap.read_count)

        logger.info("quantifying with salmon...")
        salmon = Salmon()
        expression = salmon.run(
            fasta,
            bam,
            threads=args.threads,
            output_dir=str(result_dir / "salmon"),
            error_model=not args.no_error_model,
        )

    logger.info("assigning fragments and computing read metrics...")
    read_metrics = ReadMetrics(assembly)
    read_metrics.run(
        bam,
        expression,
        fragments=snap.read_count,
        read_length=get_read_length(left),
        threads=args.threads,
        max_memory=args.max_memory,
    )
    result.update(read_metrics.read_stats())
    log_metrics("mapping metrics", result, READ_STATS_KEYS)

    optimiser = ScoreOptimiser(
        assembly=assembly,
        fragments=read_metrics.fragments,
        good=read_metrics.good,
    )
    score = optimiser.raw_score()
    weighted = optimiser.weighted_score()
    prefix = strip_gzip_suffix(assembly_path).name
    optimal, cutoff = optimiser.optimal_score(
        csv_path=str(result_dir / f"{prefix}_score_optimisation.csv")
    )
    assembly.classify_contigs(cutoff)

    result["score"] = score
    result["optimal_score"] = optimal
    result["cutoff"] = cutoff
    result["weighted"] = weighted

    logger.info("TRANSRATE ASSEMBLY SCORE     %.4f", score)
    logger.info("TRANSRATE OPTIMAL SCORE      %.4f", optimal)
    logger.info("TRANSRATE OPTIMAL CUTOFF     %.4f", cutoff)
    logger.info(
        "good contigs                 %d / %d",
        assembly.good_contigs,
        assembly.size,
    )

    write_contigs_csv(assembly, str(result_dir / CONTIGS_CSV))

    # DEFERRED_BAM_DELETE: not here. This function is one assembly of
    # possibly several, and it is followed by the run's own final step --
    # writing assemblies.csv. Deleting the BAM at the end of *this* assembly
    # meant a run that died later had already thrown away the hours that
    # produced it, and a rerun had to map again from nothing. The caller
    # deletes once the whole run has succeeded.
    if not args.keep_bam and os.path.exists(bam):
        if defer_delete is None:
            os.remove(bam)
        else:
            defer_delete.append(bam)

    # Releases the index lock; see INDEX_LOCK in pytransrate.mapper. On the
    # error paths the kernel does it when the process exits, which is the
    # next thing that happens -- an assembly that raises aborts the run.
    snap.close()

    return result


def invocation(argv=None) -> str:
    """The command line as run, quoted so it can be pasted back.

    Recorded because a log without it cannot be matched to the settings that
    produced it, which is exactly what comparing two runs requires.
    """
    parts = list(sys.argv) if argv is None else ["pytransrate", *argv]
    return " ".join(shlex.quote(str(part)) for part in parts)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.loglevel)
    if not args.no_banner:
        print_banner()
    logger.info("command: %s", invocation(argv))

    try:
        assemblies = check_arguments(args)

        output = Path(os.path.abspath(os.path.expanduser(args.output)))
        output.mkdir(parents=True, exist_ok=True)

        outfile = output / ASSEMBLIES_CSV
        if outfile.exists():
            raise CommandError(
                f"{ASSEMBLIES_CSV} would be overwritten in {output}; "
                "please choose a different output directory"
            )

        results = []
        # See DEFERRED_BAM_DELETE. Nothing in here is removed until the run
        # has produced assemblies.csv; a run that fails keeps every BAM it
        # made, so the rerun costs minutes.
        finished_bams: list[str] = []
        for assembly_path in assemblies:
            name = strip_gzip_suffix(assembly_path).stem
            result_dir = output / name if len(assemblies) > 1 else output
            result_dir.mkdir(parents=True, exist_ok=True)

            cwd = os.getcwd()
            os.chdir(result_dir)
            try:
                results.append(
                    analyse_assembly(
                        assembly_path, args, result_dir,
                        defer_delete=finished_bams,
                    )
                )
            finally:
                os.chdir(cwd)

        logger.info("writing analysis results to %s", outfile)
        write_assemblies_csv(
            results,
            str(outfile),
            with_reads=bool(args.left and args.right),
        )

        # The run has succeeded and everything derived from the BAMs is on
        # disk, so they are now genuinely spare. See DEFERRED_BAM_DELETE.
        for bam in finished_bams:
            try:
                if os.path.exists(bam):
                    logger.info("run complete; removing BAM %s", bam)
                    os.remove(bam)
                marker = align_done_path(bam)
                if marker.exists():
                    marker.unlink()
            except OSError as exc:
                logger.warning("could not remove %s (%s)", bam, exc)
    except (
        CommandError,
        AssemblyError,
        MalformedBamError,
        # plain_path, when the name the decompressed assembly wants is taken.
        FileExistsError,
    ) as exc:
        logger.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
