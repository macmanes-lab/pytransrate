"""Command-line interface.

Port of ``lib/transrate/cmdline.rb``.  The option names ORP depends on are
unchanged -- ``-a/--assembly``, ``-o/--output``, ``-t/--threads``,
``--left``, ``--right`` -- so existing invocations keep working:

    transrate -o <dir> -t <cpu> -a <fasta> --left <r1> --right <r2>

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
import sys
from pathlib import Path

from transrate import __version__
from transrate.assembly import Assembly, AssemblyError
from transrate.cmd import CommandError
from transrate.mapper import Snap
from transrate.output import write_assemblies_csv, write_contigs_csv
from transrate.quantify import Salmon
from transrate.read_metrics import ReadMetrics, get_read_length
from transrate.score import ScoreOptimiser

logger = logging.getLogger("transrate")

#: Written into the output directory, as the Ruby did.
ASSEMBLIES_CSV = "assemblies.csv"
CONTIGS_CSV = "contigs.csv"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="transrate",
        description=(
            "Analyse a de-novo transcriptome assembly using sequence-based "
            "and read-mapping metrics."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-a", "--assembly",
        required=True,
        help="assembly file(s) in FASTA format, comma-separated",
    )
    parser.add_argument("--left", help="left reads in FASTQ, comma-separated")
    parser.add_argument("--right", help="right reads in FASTQ, comma-separated")
    parser.add_argument(
        "-r", "--reference",
        help="reference proteome/transcriptome (not implemented in this port)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=8, help="number of threads"
    )
    parser.add_argument(
        "-o", "--output", default="transrate_results", help="output directory"
    )
    parser.add_argument(
        "--loglevel",
        default="info",
        choices=["error", "warn", "info", "debug"],
        help="logging verbosity",
    )
    parser.add_argument(
        "--no-error-model",
        action="store_true",
        help=(
            "don't pass --errorModel to salmon. Only sensible if your BAM "
            "carries AS tags; snap-aligner does not emit them"
        ),
    )
    parser.add_argument(
        "--keep-bam",
        action="store_true",
        help="keep the alignment BAM instead of deleting it on success",
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)5s] %(asctime)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


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


def analyse_assembly(assembly_path, args, result_dir: Path) -> dict:
    """Run every requested metric for one assembly."""
    logger.info("loading assembly: %s", assembly_path)
    assembly = Assembly(assembly_path)

    result = {"assembly": str(assembly_path)}
    logger.info("calculating contig metrics...")
    result.update(assembly.basic_stats())
    result.update(assembly.contig_metrics())

    with_reads = bool(args.left and args.right)
    if not with_reads:
        logger.info("no reads provided, skipping read diagnostics")
        write_contigs_csv(
            assembly, str(result_dir / CONTIGS_CSV), with_reads=False
        )
        return result

    left, right = args.left, args.right  # already absolute, see check_arguments

    logger.info("mapping reads with snap-aligner...")
    snap = Snap()
    snap.build_index(assembly_path, threads=args.threads)
    bam = snap.map_reads(left, right, threads=args.threads)
    logger.info("%d fragments in library", snap.read_count)

    logger.info("quantifying with salmon...")
    salmon = Salmon()
    expression = salmon.run(
        assembly_path,
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
    )
    result.update(read_metrics.read_stats())

    optimiser = ScoreOptimiser(
        assembly=assembly,
        fragments=read_metrics.fragments,
        good=read_metrics.good,
    )
    score = optimiser.raw_score()
    weighted = optimiser.weighted_score()
    prefix = Path(assembly_path).name
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

    if not args.keep_bam and os.path.exists(bam):
        os.remove(bam)

    return result


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.loglevel)

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
        for assembly_path in assemblies:
            name = Path(assembly_path).stem
            result_dir = output / name if len(assemblies) > 1 else output
            result_dir.mkdir(parents=True, exist_ok=True)

            cwd = os.getcwd()
            os.chdir(result_dir)
            try:
                results.append(
                    analyse_assembly(assembly_path, args, result_dir)
                )
            finally:
                os.chdir(cwd)

        logger.info("writing analysis results to %s", outfile)
        write_assemblies_csv(
            results,
            str(outfile),
            with_reads=bool(args.left and args.right),
        )
    except (CommandError, AssemblyError) as exc:
        logger.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
