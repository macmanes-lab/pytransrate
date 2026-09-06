"""Read mapping with snap-aligner.

Port of ``lib/transrate/snap.rb``, targeting **snap-aligner 2.0.5**
(amplab/snap) rather than the abandoned Blahah 1.0dev.96 fork.

Every flag the Ruby passed still exists in 2.0.5 and was verified against the
real binary: index takes ``-s``, ``-t``, ``-bSpace`` and ``-locationSize``;
paired takes ``-s min max``, ``-H``, ``-h``, ``-d``, ``-t``, ``-b``, ``-M``,
``-D``, ``-om``, ``-omax`` and ``-mcp``.  So the command line is unchanged.

What is *not* unchanged is the output: see SOFT_CLIPPING below.
"""

from __future__ import annotations

import logging
import os
import shutil
import re
from pathlib import Path

from transrate.cmd import CommandError, run, which

__all__ = ["SnapError", "Snap"]

logger = logging.getLogger("transrate")


class SnapError(CommandError):
    """snap-aligner failed."""


# ---------------------------------------------------------------------------
# SOFT_CLIPPING
#
# snap-aligner 2.0.0 introduced soft clipping, and 2.0.5 applies it by
# default ("when a read (or pair) doesn't align, try soft clipping the read
# (or pair) to find an alignment.  This is the default.").  Measured on reads
# carrying a 25bp foreign prefix, 180 of 407 alignments came back clipped
# (25S75M, 24S76M, 23S2M1D75M, ...).
#
# The old bam-read binary advanced its reference cursor on 'S' as though it
# were a deletion, so those alignments would have had their coverage
# displaced rightward by the clip length.  transrate.bam_metrics follows the
# SAM spec instead (see SOFT_CLIP_FIX there).  Do not reintroduce bam-read
# alongside this aligner.
#
# Note also that snap-aligner emits **no AS tag** -- only NM -- which
# constrains both salmon's scoring and our own fragment assignment.  See
# transrate.quantify and transrate.assign.
# ---------------------------------------------------------------------------

#: locationSize values tried when the index overflows, as in the Ruby.
#: SNAP enforces 4..8 inclusive and hard-exits outside it.
_LOCATION_SIZES = range(4, 9)

#: Written by snap-aligner once an index is complete. The directory alone is
#: not proof of a usable index -- a build that died partway leaves the
#: directory behind with some of its files.
_INDEX_MARKER = "GenomeIndex"

_OVERFLOW_PATTERNS = (
    re.compile(r"Ran out of overflow table namespace"),
    re.compile(r"Trying to use too many overflow entries"),
)

_UNMATCHED_IDS = re.compile(r"Unmatched\s+read\s+IDs\s+(.*?)\s+and\s+(.*?)Use", re.S)


class Snap:
    """Builds a snap index and maps paired reads against it."""

    def __init__(self, binary: str | None = None):
        self.binary = binary or which("snap-aligner")
        self.index_name: str | None = None
        self.index_built = False
        self.bam: str | None = None
        self.read_count: int = 0
        self._read_count_file: str | None = None

    # -- index ------------------------------------------------------------

    def build_index(
        self,
        fasta,
        threads: int = 8,
        seed_size: int = 23,
        location_size: int | None = None,
    ) -> str:
        """Build a snap index.

        ``-locationSize`` sets how many bytes each genome location occupies.
        Four is enough for most assemblies, but a large or highly repetitive
        one exhausts the overflow table and snap refuses to index it:

            Ran out of overflow table namespace. This genome cannot be
            indexed with this seed and location size.  Increase at least one.

        Args:
            fasta: assembly to index.
            threads: passed as ``-t`` (no space, as snap requires).
            seed_size: passed as ``-s``.
            location_size: fix ``-locationSize`` at this value and do not
                sweep. ``None`` (the default) tries 4 and steps up to 8 on
                an overflow error, as the Ruby did. Fixing it is worth doing
                when you already know an assembly needs a larger value --
                each failed attempt is a full index build.

        Raises:
            SnapError: if the build fails, or overflows at every size tried.
        """
        fasta = Path(fasta)
        self.index_name = fasta.stem
        index_dir = Path(self.index_name)

        # Require the marker, not just the directory: a build killed partway
        # leaves the directory behind, and trusting it yields a corrupt index.
        if (index_dir / _INDEX_MARKER).exists():
            self.index_built = True
            return self.index_name

        if location_size is not None:
            if location_size not in _LOCATION_SIZES:
                raise SnapError(
                    f"location_size must be between {_LOCATION_SIZES.start} and "
                    f"{_LOCATION_SIZES.stop - 1} inclusive, got {location_size}"
                )
            sizes = [location_size]
        else:
            sizes = list(_LOCATION_SIZES)

        last_error = ""
        for size in sizes:
            args = [
                self.binary, "index", str(fasta), self.index_name,
                "-s", seed_size,
                f"-t{threads}",
                "-bSpace",              # contig name ends at the first space
                "-locationSize", size,
            ]
            result = run(args)
            if result.ok:
                self.index_built = True
                return self.index_name

            last_error = result.stderr or result.stdout
            if any(p.search(last_error) for p in _OVERFLOW_PATTERNS):
                # Clear the partial index before retrying. rmdir() will not
                # do -- snap leaves Genome, GenomeIndex, GenomeIndexHash and
                # OverflowTable behind, and rmdir only removes empty dirs.
                if index_dir.is_dir():
                    shutil.rmtree(index_dir, ignore_errors=True)
                if size == sizes[-1]:
                    break
                logger.warning(
                    "snap index overflowed at -locationSize %d, retrying at %d",
                    size,
                    size + 1,
                )
                continue
            raise SnapError(f"Failed to build snap index\n{last_error}")

        hint = (
            " Every location size from "
            f"{sizes[0]} to {sizes[-1]} overflowed; try a larger --seed-size."
            if len(sizes) > 1
            else f" Retry without --location-size to sweep {_LOCATION_SIZES.start}"
                 f"-{_LOCATION_SIZES.stop - 1}, or use a larger --seed-size."
        )
        raise SnapError(f"Failed to build snap index.{hint}\n{last_error}")

    # -- mapping ----------------------------------------------------------

    def build_paired_command(self, left, right, threads: int, output: str):
        """Assemble the ``snap-aligner paired`` command.

        Kept identical to the Ruby's flags; all were verified present in
        snap-aligner 2.0.5.
        """
        args = [self.binary, "paired", self.index_name]
        for l, r in zip(str(left).split(","), str(right).split(",")):
            args += [l, r]
        args += [
            "-o", output,
            "-s", 0, 1000,        # min/max spacing between paired-read starts
            "-H", 300000,         # max seed hits in paired mode
            "-h", 2000,           # max seed hits when reverting to single
            "-d", 30,             # max edit distance
            "-t", threads,
            "-b",                 # bind threads to cores
            "-M",                 # M-style CIGAR (now the default, kept explicit)
            "-D", 5,              # extra search depth, needed for -om
            "-om", 5,             # report multiple alignments
            "-omax", 10,          # cap alignments per pair
            "-mcp", 10000000000000,
        ]
        return args

    def map_reads(self, left, right, threads: int = 8, output=None) -> str:
        """Map paired reads, returning the path to the BAM.

        The BAM is left in **read order**, not coordinate sorted: both
        :func:`transrate.assign.assign_fragments` and
        :func:`transrate.bam_metrics.estimate_realistic_distance` require
        mates and multi-mappings to be adjacent.
        """
        if not self.index_built:
            raise SnapError("Index not built")

        lbase = Path(str(left).split(",")[0]).name
        rbase = Path(str(right).split(",")[0]).name
        index = Path(self.index_name).name
        self.bam = str(Path(output or f"{lbase}.{rbase}.{index}.bam").resolve())
        self._read_count_file = f"{lbase}-{rbase}-read_count.txt"

        if os.path.exists(self.bam):
            self._load_read_count(left)
            return self.bam

        args = self.build_paired_command(left, right, threads, self.bam)
        result = run(args)
        self._save_read_count(result.stdout)
        self._save_logs(result.stdout, result.stderr)

        if not result.ok:
            match = _UNMATCHED_IDS.search(result.stderr)
            if match:
                raise SnapError(
                    "snap found unmatched read IDs in the input fastq files.\n"
                    f"Left files contained read id\n{match.group(1).strip()}\n"
                    f"and right files contained read id\n{match.group(2).strip()}\n"
                    "at the same position in the file."
                )
            raise SnapError(f"snap failed\n{result.stderr}")

        return self.bam

    # -- read counting ----------------------------------------------------

    def _save_logs(self, stdout: str, stderr: str) -> None:
        Path("logs").mkdir(exist_ok=True)
        with open("logs/snap.log", "a") as handle:
            handle.write(stdout)
            handle.write(stderr)

    def _save_read_count(self, stdout: str) -> None:
        """Pull the total read count out of snap's summary table.

        snap prints a header row then a row of figures whose first column is
        the total reads processed, comma-grouped. Fragments are half that.
        """
        for line in stdout.splitlines():
            cols = line.split()
            if len(cols) > 5 and re.fullmatch(r"[0-9,]+", cols[0]):
                self.read_count = int(cols[0].replace(",", "")) // 2
                if self._read_count_file:
                    with open(self._read_count_file, "w") as handle:
                        handle.write(f"{self.read_count}\n")

    def _load_read_count(self, reads) -> None:
        """Recover the read count when reusing an existing BAM."""
        self.read_count = 0
        if self._read_count_file and os.path.exists(self._read_count_file):
            with open(self._read_count_file) as handle:
                self.read_count = int(handle.read().strip() or 0)
            return

        for path in str(reads).split(","):
            try:
                with open(path, "rb") as handle:
                    lines = sum(1 for _ in handle)
                self.read_count += lines // 4
            except OSError:
                logger.warning("couldn't count reads in %s", path)

        if self._read_count_file:
            with open(self._read_count_file, "w") as handle:
                handle.write(f"{self.read_count}\n")
