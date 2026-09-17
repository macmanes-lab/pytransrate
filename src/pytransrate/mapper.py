"""Read mapping with snap-aligner.

Port of ``lib/transrate/snap.rb``, targeting **snap-aligner 2.0.5**
(amplab/snap) rather than the abandoned Blahah 1.0dev.96 fork.

Every flag the Ruby passed still exists in 2.0.5 and was verified against the
real binary: index takes ``-s``, ``-t``, ``-bSpace``, ``-p`` and
``-locationSize``; paired takes ``-s min max``, ``-H``, ``-h``, ``-d``,
``-t``, ``-b``, ``-M``, ``-D``, ``-om`` and ``-omax``.  ``-mcp`` is no longer passed by default --
see MAX_CANDIDATE_POOL below.

What is *not* unchanged is the output: see SOFT_CLIPPING below.
"""

from __future__ import annotations

import datetime
import fcntl
import logging
import os
import shutil
import re
import signal
from pathlib import Path

from pytransrate.cmd import CommandError, run, which
from pytransrate.compression import open_binary

__all__ = ["SnapError", "Snap", "bam_is_complete", "alignment_is_done",
           "align_done_path"]

logger = logging.getLogger("pytransrate")


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
# displaced rightward by the clip length.  pytransrate.bam_metrics follows the
# SAM spec instead (see SOFT_CLIP_FIX there).  Do not reintroduce bam-read
# alongside this aligner.
#
# Note also that snap-aligner emits **no AS tag** -- only NM -- which
# constrains both salmon's scoring and our own fragment assignment.  See
# pytransrate.quantify and pytransrate.assign.
# ---------------------------------------------------------------------------

#: locationSize values tried when the index overflows, as in the Ruby.
#: SNAP enforces 4..8 inclusive and hard-exits outside it.
_LOCATION_SIZES = range(4, 9)

#: Written by snap-aligner once an index is complete. The directory alone is
#: not proof of a usable index -- a build that died partway leaves the
#: directory behind with some of its files.
_INDEX_MARKER = "GenomeIndex"

#: Everything snap-aligner says, indexing and mapping alike, as it says it.
#: Relative: the CLI chdirs into the result directory, so each assembly gets
#: its own.
LOG_PATH = Path("logs") / "snap.log"

# ---------------------------------------------------------------------------
# INDEX_LOCK
#
# The index lives in a directory named after the assembly, so two runs of the
# same assembly into the same output directory share it.  Only one code path
# ever deletes an index -- build_index clearing a partial build before
# retrying at a larger -locationSize -- but that is enough: if the second run
# reaches it while the first is already mapping, snap-aligner loses its
# genome mid-alignment.
#
# So an exclusive lock is taken before the directory is touched and held
# until the Snap object is closed, which covers mapping as well as building.
# flock is released by the kernel when the holder dies, so a killed run
# cannot leave a lock behind for the next one to trip over -- which a lock
# file tested with O_EXCL would.
#
# The lock sits beside the index rather than inside it, because rmtree would
# otherwise delete the file out from under its own holder and the next
# process would happily lock a fresh inode.
# ---------------------------------------------------------------------------

_LOCK_SUFFIX = ".index.lock"

# ---------------------------------------------------------------------------
# MAX_CANDIDATE_POOL
#
# snap.rb passed `-mcp 10000000000000` with the comment "increase mcp to silly
# high value to dec incidence of common SNAP fail" -- a workaround for
# snap-aligner 1.0dev.96, the abandoned fork this port no longer targets.
#
# It never meant what it looks like. snap parses it with atoi() into an int,
# and 10^13 is far past INT_MAX (2147483647), so the conversion is undefined
# behaviour. Measured: atoi("10000000000000") returns 1316134912 -- a wrapped
# garbage value that varies by platform, libc and optimisation level. The
# request for "effectively unlimited" silently became an arbitrary number.
#
# snap has a sane DEFAULT_MAX_CANDIDATE_POOL_SIZE, so -mcp is no longer passed
# unless asked for. Pass a value explicitly to restore it, and keep it under
# INT_MAX if you do.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# MULTI_ALIGNMENT_SETTINGS
#
# snap-aligner 2.0.5 dies with SIGFPE on a real ORP assembly under the flags
# snap.rb used (-H 300000 -D 5 -om 5 -omax 10). Bisected on the failing data:
# removing -mcp did not help, nor did dropping -H to 4000; removing
# -D/-om/-omax did. So the fault is in the multiple-alignment path, reached
# only at a scale and redundancy this project's synthetic tests do not hit --
# 27k redundant contigs with 116k secondary alignments still aligns cleanly
# here.
#
# That diagnosis was half right, and SNAP_171 below has the rest of it: the
# fault is indeed in the multiple-alignment path, but the flags are not what
# provokes it, and the defaults below do NOT make a large enough run safe.
# They remain the right settings on their own merits, measured below.
#
# The defaults below are the configuration that runs on that data. Two of
# them differ from the Ruby:
#
#   -om 5 -> 2   "accept alignments within edit distance 5 of the best" is far
#                wider than a redundant transcriptome needs, where duplicate
#                contigs differ by 0-2 bases. Measured on a redundant 27k-contig
#                assembly, tightening 5 -> 1 retained 98.6% of secondary
#                alignments and 99.98% of fragments hitting >1 contig, so the
#                multi-mapping signal the assignment step needs survives intact.
#   -D 5 -> 2    extra search depth; must track -om.
#
# -mpc 1 is new, and is the right shape for this job independently of the
# crash: it caps alignments per contig, applied before -omax, so a fragment
# contributes its best placement on each candidate contig rather than several
# placements on one. That is what transrate.assign consumes, and it is the
# precondition accumulate_metrics states in its docstring.
#
# MEASURED, after the fact, on SRR1789336 (28,976,658 fragments) against
# three assemblies at -mpc 1, 2 and 0. Loosening the cap raises the transrate
# score every time, and the rise is an artefact of counting rather than a
# better assembly:
#
#                             score    from good-rate   from contig geomean
#    ORP  (101,342 contigs)  +0.0050   +0.0030 ( 60%)   +0.0020 ( 40%)
#    s55  (124,787 contigs)  +0.0027   +0.0030 (112%)   -0.0003 (-12%)
#    s75  ( 99,050 contigs)  +0.0031   +0.0036 (116%)   -0.0005 (-16%)
#
# (score = geomean(contig scores) * good_mappings/fragments; see score.py.
# -mpc 0 shown. -mpc 2 lands within 0.0005 of it on all three, so raising the
# cap to 2 buys nothing. That gap is not quite zero -- -mpc 2 scores 0.0004
# above -mpc 0 on ORP, same sign across both replicates, ~6x the noise floor
# below -- but it is an order of magnitude smaller than the jump off -mpc 1,
# so the cap behaves almost binarily at 1 vs >1.)
#
# The effect is real. The whole ORP matrix was run twice, giving a replicate
# at every setting: score reproduces to 2e-5 (8e-5 at -mpc 2), optimal_score
# to 1e-5, fragments_mapped exactly, and good_mappings to within 1,800 of
# 24.5M. So the +0.0050 above is 60-250x the run-to-run noise, and the
# per-contig counts quoted below reproduce to ~0.05% (sCnuc 7749/19690 on the
# first pass, 7747/19684 on the second). snap's documented nondeterminism
# (amplab/snap#72) is real but far too small to matter here.
#
# But the effect points the wrong way. -mpc caps alignments per contig, so it
# cannot make a previously unalignable fragment align: the set of fragments
# carrying an alignment is identical across the three runs. fragments_mapped
# nonetheless rose by 205k-277k, and since it increments once per read-1
# record on the assigned contig, the only thing that can move it is one
# fragment being counted several times.
# good_mappings and bad_mappings rise together, where a reassignment would
# trade one for the other, and good_mappings/fragments -- an inflated
# numerator over a fixed denominator -- carries most or all of the score
# change.
#
# The two terms behave nothing alike, which is the tell. The good-rate term
# is ~+0.003 on all three assemblies, spanning 99k-125k contigs and a 6x
# range of uncovered bases, and fragments_mapped moves by 0.7-1.0% of the
# library every time: a property of the reads, not of what they were mapped
# to. The geomean term swings with the assembly instead, and does so in the
# order the mechanism predicts -- sCcov can only gain where there was
# uncovered sequence to reclaim:
#
#    ORP  10.8% bases uncovered  ->  sCcov moved on 6.6% of contigs, +0.0020
#    s55   3.3%                  ->                  1.7%,           -0.0003
#    s75   1.7%                  ->                  1.2%,           -0.0005
#
# So ORP is the exception rather than the pattern: it had enough uncovered
# sequence for the inflated coverage to lift the geomean. On a normally
# covered assembly the contigs measurably worsen and the headline score rises
# anyway. sCnuc falls decisively everywhere (s55: 2,638 contigs up against
# 16,789 down; s75: 1,832 against 16,909), because the extra placements are
# worse than the ones already there.
#
# A second-order harm comes with it. With -omax 10 fixed, -mpc 1 gives a
# fragment up to ten alignments across ten distinct contigs; -mpc 0 lets all
# ten land on one. Loosening the cap narrows the candidate set the assignment
# step chooses from rather than widening it.
#
# So do not raise this default. Note also that it is not a workaround for the
# records htslib refuses -- those appear at every setting (0, 1 and 3 skipped
# at -mpc 1, 2 and 0 on the ORP assembly); see MALFORMED_RECORDS in
# pytransrate.bam_metrics. scripts/compare_transrate_runs.py reproduces the
# comparison above from any set of run directories.
#
# -H drops to snap's own default of 4000. The Ruby's 300000 was 75x that,
# with no rationale recorded, and it sizes the scoring candidate pool
# (scoringCandidatePoolSize = min(mcp, maxBigHits * maxSeeds * 2)).
#
# All are overridable; see the --max-seed-hits family in pytransrate.cli.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# SNAP_171
#
# snap 2.0.x dies with SIGFPE partway through alignment on large runs, and
# no setting this port exposes prevents it: amplab/snap#171.
#
# A read reaches computeGlobalScore with patternLen == 0, so numVec == 0
# (AffineGapVectorized.cpp:186) and line 351 divides by it. Upstream traced
# the zero in July 2025: when the writer fills its output buffer partway
# through writing a read's alignments it flushes, takes a fresh buffer and
# retries the write, and the back-clipping from a secondary alignment was
# wrongly retained across the retry. A read whose secondary alignment clipped
# it to exactly half its length is then clipped to nothing on the retry.
# Any other length gave a silently wrong alignment instead of a crash.
#
# So the trigger is secondary alignments (-om/-omax) plus a buffer refill at
# the wrong instant, which makes it a function of OUTPUT VOLUME, not of
# assembly size -- it merely looks like a large-assembly bug because that is
# what produces a large BAM. On a ~160 GB BAM the writer refills more or less
# continuously and a per-read-rare coincidence becomes a certainty: measured
# on a 5,354,958-contig, 5.66 Gbp merged assembly, snap dies ~17 minutes in
# with the BAM past 100 GB, on every attempt. The same flags on smaller
# assemblies run clean, which is why no test here catches it.
#
# Tuning does not help. The original report used -H 300000 -D 5 -om 5; this
# port uses -H 4000 -D 2 -om 2 and still dies. Only dropping -om/-omax
# entirely avoids it, by removing the secondary alignments the bug needs --
# which is not available to us, since those are precisely what
# pytransrate.assign consumes.
#
# Fixed upstream in commit 0e0997b, released as 2.0.6.dev.2 on snap's dev
# branch. It is not in master and not in the bioconda 2.0.5 package, so a
# large assembly needs a dev build. _how_it_died says so on SIGFPE.
# ---------------------------------------------------------------------------

#: snap's -H, max hits for the intersecting aligner. snap's own default.
DEFAULT_MAX_SEED_HITS = 4000

#: snap's -d, maximum edit distance per read or pair.
DEFAULT_EDIT_DISTANCE = 30

#: snap's -D, extra search depth. Must be >= the -om value.
DEFAULT_EXTRA_SEARCH_DEPTH = 2

#: snap's -om, extra edit distance admitted for secondary alignments.
DEFAULT_MULTI_EDIT_DISTANCE = 2

#: snap's -omax, cap on secondary alignments per pair.
DEFAULT_MAX_ALIGNMENTS_PER_PAIR = 10

#: snap's -mpc, cap on alignments per contig, applied before -omax.
DEFAULT_MAX_ALIGNMENTS_PER_CONTIG = 1

# ---------------------------------------------------------------------------
# CONTIG_PADDING
#
# snap pads every contig with Ns so an alignment cannot run off one contig
# into the next, and those Ns are real bases as far as the index is
# concerned: FASTA.cpp sizes the genome as fileSize + (nContigs + 1) *
# padding and writes the padding in, so getCountOfBases() counts it and the
# -locationSize check below is applied to the total.
#
# snap's own default is 2000 (GenomeIndex.cpp DEFAULT_PADDING), chosen for
# genomes, where a few hundred contigs make the padding a rounding error.  A
# transcriptome inverts that: the size the location namespace has to cover is
#
#     n_bases + padding * (n_contigs + 1)
#
# and a de-novo assembly has contigs by the million, averaging a kilobase or
# two.  At 2000 the padding can exceed the assembly -- 1.5M contigs contribute
# 3 Gbp of Ns -- which is enough on its own to cross the 4-byte location
# ceiling of 2**32 - 16 bases and force the sweep below up to -locationSize 5.
# Every genome location then costs 5 bytes rather than 4, for an index the
# aligner has to hold in memory for the length of the run.  On the three of
# the four messages below that snap only reaches after building, the failed
# attempt costs a full index build as well; on "Genome is too big" it does
# not, since that one is checked up front.
#
# 1000 is the largest reduction that costs nothing on either count snap
# documents for this value:
#
#   "This must be as large as the largest edit distance you'll ever use, and
#    there's a performance advantage to have it be bigger than any read you'll
#    process or gap between paired-end reads."
#
# The edit distance is a correctness floor and is DEFAULT_EDIT_DISTANCE, 30 --
# two orders of magnitude clear.  The second clause is a performance note, not
# a correctness one, and the gap it refers to is the -s maximum in
# build_paired_command, which is 1000.  Padding of 1000 sits exactly at that
# bound rather than above it, so a pair whose ends straddle two adjacent
# contigs is no longer separated by more than the maximum spacing.  It cannot
# be called a proper pair regardless -- crossing the padding means crossing
# 1000 Ns, which no alignment within an edit distance of 30 survives -- so
# what is at stake is snap doing the work to reject it, not the rejection.
# Raise this above 1000 + read length if that ever shows up in a profile.
#
# This does change alignments, and so scores, on any assembly where it changes
# the index: contigs sit at different genome locations.  See the note in
# CHANGELOG.md.
# ---------------------------------------------------------------------------

#: snap index -p, Ns inserted between contigs. Below snap's own default of
#: 2000; see CONTIG_PADDING.
DEFAULT_PADDING = 1000

# ---------------------------------------------------------------------------
# LOCATION_SIZE_FAILURES
#
# snap-aligner refuses to index a genome that does not fit the location-size
# namespace in four different places, and all four are cured by a larger
# -locationSize. Matching only some of them means the sweep below never runs
# and the user is told to do by hand what pytransrate is meant to do for them,
# so all four are listed here (GenomeIndex.cpp, snap 2.0.5):
#
#   "Genome is too big for %d byte genome locations.  Specify a larger
#    location size with -locationSize"
#       countOfBases > 2**(locationSize*8) - 16, checked before any index
#       work -- the cheapest of the four to hit, and the one a large merged
#       assembly hits first.
#
#   "Ran out of overflow table namespace. This genome cannot be indexed with
#    this seed and location size.  Increase at least one."
#
#   "Trying to use too many overflow entries.  To index this genome, you
#    either need a larger seed size or a larger location size."
#
#   "Not enough address space to index this genome with this seed size.  Try
#    a larger seed or location size."
#       despite the wording this is not a RAM limit: the bound is
#       InvalidGenomeLocation, which is 2**(locationSize*8) - 1, so a larger
#       location size does raise it.
#
# Match the distinctive head of each message only. snap wraps and punctuates
# these inconsistently (note the double spaces), and the trailing advice is
# what upstream is most likely to reword.
# ---------------------------------------------------------------------------

_LOCATION_SIZE_PATTERNS = (
    re.compile(r"Genome is too big for \d+ byte genome locations"),
    re.compile(r"Ran out of overflow table namespace"),
    re.compile(r"Trying to use too many overflow entries"),
    re.compile(r"Not enough address space to index this genome"),
)

_UNMATCHED_IDS = re.compile(r"Unmatched\s+read\s+IDs\s+(.*?)\s+and\s+(.*?)Use", re.S)

# ---------------------------------------------------------------------------
# SILENT_OPTION_REJECTION
#
# snap-aligner prints "Didn't understand options starting at ..." followed by
# its usage text, writes no BAM -- and exits 0. A caller that trusts the exit
# code sees success and carries on with a path to a file that does not exist,
# failing later and somewhere unrelated. So the mapping step checks the
# output itself rather than the return code alone.
# ---------------------------------------------------------------------------

_BAD_OPTIONS = re.compile(r"Didn't understand options starting at (.*)")


#: snap's ceiling on a 4-byte location: InvalidGenomeLocation, less the
#: sentinel slack GenomeIndex.cpp leaves above it.
_FOUR_BYTE_CEILING = 2 ** 32 - 16


def _count_contigs(fasta: Path) -> int:
    """'>' at the start of a line, counted without holding the file."""
    count = 0
    last = b"\n"
    with open(fasta, "rb") as handle:
        while chunk := handle.read(1 << 22):
            if last == b"\n" and chunk[:1] == b">":
                count += 1
            count += chunk.count(b"\n>")
            last = chunk[-1:]
    return count


def _padding_report(fasta: Path, padding: int) -> str:
    """How much of the genome snap is measuring is padding, not assembly.

    snap sizes the genome as ``fileSize + (nContigs + 1) * padding``
    (FASTA.cpp) and applies the location-size ceiling to that total, so on a
    fragmented transcriptome the padding is routinely several times the
    assembly.  Without this the user sees only that four bytes were not
    enough, and reaches for --seed-size or a bigger machine when --padding is
    the term actually driving the number.

    Computed only on the failure path: it reads the whole FASTA, which is
    nothing beside the index build it is explaining, but is not free.
    """
    try:
        file_size = fasta.stat().st_size
        contigs = _count_contigs(fasta)
    except OSError:
        return ""

    pad_bases = (contigs + 1) * padding
    genome = file_size + pad_bases
    if not genome:
        return ""
    return (
        f"snap measures this genome as {file_size / 1e9:.2f} GB of FASTA "
        f"(it counts the file's bytes, headers and newlines included) plus "
        f"{contigs:,} contigs x {padding} bp of padding "
        f"({pad_bases / 1e9:.2f} Gbp) = {genome / 1e9:.2f} Gbp, "
        f"{100.0 * pad_bases / genome:.0f}% of it padding; the 4-byte "
        f"ceiling is {_FOUR_BYTE_CEILING / 1e9:.2f} Gbp. Padding is inert Ns "
        f"and still costs index and memory, so --padding is usually the "
        f"cheaper lever: its floor is the read length plus --edit-distance, "
        f"below which an alignment could cross from one contig into the next."
    )


def _how_it_died(returncode: int) -> str:
    """Describe a return code, naming the signal when there was one.

    subprocess reports a signal death as a negative return code, and "-8" on
    its own tells a user nothing.  The distinction is worth spelling out:
    a non-zero exit is snap declining to do something and usually saying why,
    whereas a signal is snap crashing, where the flags that provoked it are
    the only evidence there is.
    """
    if returncode >= 0:
        return f"exit {returncode}"
    try:
        name = signal.Signals(-returncode).name
    except ValueError:
        name = "unrecognised signal"
    note = ""
    if -returncode == signal.SIGFPE:
        # SNAP_171. Not "usually" this: snap 2.0.x has exactly one known
        # divide-by-zero, and nothing this port exposes avoids it, so
        # sending the user round a tuning loop wastes hours of their time.
        note = (
            " -- this is amplab/snap#171 until proved otherwise: a "
            "divide-by-zero in snap 2.0.x reached when the writer refills "
            "its output buffer mid-read, so it tracks the size of the BAM "
            "rather than the size of the assembly. No pytransrate setting "
            "avoids it, and lowering the multi-alignment flags does not "
            "either. Fixed upstream in 2.0.6.dev.2, which must be built "
            "from snap's dev branch. See SNAP_171 in mapper.py"
        )
    elif -returncode in (signal.SIGSEGV, signal.SIGBUS, signal.SIGILL):
        note = (
            " -- snap crashed rather than reporting an error, so the last "
            "lines of the log are the evidence: they say whether it got as "
            "far as loading the index or died partway through the reads. On "
            "a large assembly this is usually the multiple-alignment path, "
            "so try lowering --max-alignments-per-pair and --max-seed-hits"
        )
    elif -returncode in (signal.SIGKILL, signal.SIGTERM):
        note = (
            " -- something outside snap stopped it; check the job's memory "
            "ceiling and wall clock rather than its flags"
        )
    return f"killed by signal {-returncode} ({name}){note}"


def _tail(text: str, lines: int = 40) -> str:
    """The last few lines of a command's output.

    snap prints a progress table and a version banner before it says what
    went wrong, and an exception carrying all of it buries the message. The
    whole thing is in LOG_PATH either way.
    """
    kept = text.strip().splitlines()[-lines:]
    return "\n".join(kept)


#: The 28-byte empty BGZF block that ends every complete BAM (SAM spec 4.1).
#: A BAM whose last bytes are not this one was still being written when
#: whatever was writing it stopped.
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003" + "00" * 9)


def bam_is_complete(path) -> bool:
    """Whether ``path`` is a BAM that was written all the way to the end.

    Resuming a run reuses the BAM it finds, and mapping a real library takes
    hours, so this is the check that decides between hours saved and metrics
    computed off half a library.  A BAM left by a killed run is the common
    case here, not a corner: this module already exists partly because snap
    gets killed by OOM killers and wall clocks.
    """
    try:
        if os.path.getsize(path) < len(_BGZF_EOF):
            return False
        with open(path, "rb") as handle:
            handle.seek(-len(_BGZF_EOF), os.SEEK_END)
            return handle.read(len(_BGZF_EOF)) == _BGZF_EOF
    except OSError:
        return False


# ---------------------------------------------------------------------------
# ALIGN_DONE
#
# Mapping a real library takes hours and the BAM it produces can be hundreds
# of gigabytes, so every decision about that file is a decision about whether
# somebody's afternoon survives. Two things used to put it at risk.
#
# The first is inference. bam_is_complete reads the BGZF end-of-file marker,
# which says the file was closed cleanly -- but not which command wrote it,
# nor whether the reads or the index have changed since. A marker file written
# only after snap has exited 0 and passed every check in map_reads says that
# directly, and records the command so a reuse can be matched to the settings
# that produced it. The EOF check is still made: the marker says the run
# finished, the marker plus the EOF says the file did too.
#
# The second is overwriting. A crashed snap leaves a large partial BAM, and
# re-mapping used to write straight over it. That file cannot be used for
# metrics -- half a library gives half the coverage -- but it is the only
# evidence of what the aligner did before it died, which is exactly what is
# wanted when the crash is amplab/snap#171 (see SNAP_171). It is now moved
# aside rather than destroyed. At most one is kept, so the cost is bounded at
# one extra BAM rather than growing with every retry.
# ---------------------------------------------------------------------------

#: Written beside the BAM once mapping has finished and been checked.
ALIGN_DONE_SUFFIX = ".align.done"

#: Where a partial BAM is moved before re-mapping, rather than overwritten.
PARTIAL_SUFFIX = ".partial"


def align_done_path(bam) -> Path:
    """The completion marker that belongs to ``bam``."""
    return Path(str(bam) + ALIGN_DONE_SUFFIX)


def alignment_is_done(bam) -> bool:
    """Whether ``bam`` was written by a mapping run that finished.

    Both halves matter: the marker says snap returned and passed its checks,
    the BGZF end-of-file marker says the file itself was closed. A BAM with
    neither is a killed run; a BAM with the EOF but no marker predates the
    marker and is reused on the strength of the EOF alone.
    """
    return align_done_path(bam).exists() and bam_is_complete(bam)


class Snap:
    """Builds a snap index and maps paired reads against it.

    Holds an exclusive lock on the index directory from :meth:`build_index`
    until :meth:`close`; see INDEX_LOCK. Usable as a context manager.
    """

    def __init__(self, binary: str | None = None):
        self.binary = binary or which("snap-aligner")
        self.index_name: str | None = None
        self.index_built = False
        self.bam: str | None = None
        self.read_count: int = 0
        self._read_count_file: str | None = None
        self._lock = None

    def close(self) -> None:
        """Release the index lock. See INDEX_LOCK."""
        if self._lock is not None:
            self._lock.close()   # closing drops the flock
            self._lock = None

    def __enter__(self) -> "Snap":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # -- index ------------------------------------------------------------

    def _lock_index(self, index_dir: Path) -> None:
        """Claim the index directory for this process. See INDEX_LOCK."""
        if self._lock is not None:
            return
        lock_path = index_dir.parent / f"{index_dir.name}{_LOCK_SUFFIX}"
        # Opened without truncating, so the holder's pid is still readable
        # when the lock is refused.
        handle = os.fdopen(os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o644), "r+")
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            held_by = handle.read().strip() or "another process"
            handle.close()
            raise SnapError(
                f"the snap index {index_dir} is in use by pid {held_by} "
                f"(lock: {lock_path}). Two runs of the same assembly in one "
                f"output directory share an index, and the second can delete "
                f"it while the first is still mapping. Give this run its own "
                f"-o directory, or wait for the other to finish."
            ) from None
        handle.seek(0)
        handle.truncate()
        handle.write(f"{os.getpid()}\n")
        handle.flush()
        self._lock = handle

    def build_index(
        self,
        fasta,
        threads: int = 8,
        seed_size: int = 23,
        location_size: int | None = None,
        padding: int = DEFAULT_PADDING,
    ) -> str:
        """Build a snap index.

        ``-locationSize`` sets how many bytes each genome location occupies.
        Four is enough for most assemblies, but a large or highly repetitive
        one exhausts that namespace and snap refuses to index it:

            Genome is too big for 4 byte genome locations.  Specify a
            larger location size with -locationSize

        See LOCATION_SIZE_FAILURES above for the four ways snap says this.

        Args:
            fasta: assembly to index.
            threads: passed as ``-t`` (no space, as snap requires).
            seed_size: passed as ``-s``.
            location_size: fix ``-locationSize`` at this value and do not
                sweep. ``None`` (the default) tries 4 and steps up to 8 when
                snap runs out of location namespace, as the Ruby did. Fixing
                it is worth doing when you already know an assembly needs a
                larger value -- each failed attempt is a full index build.
            padding: passed as ``-p`` (no space, as snap requires).  Ns
                inserted between contigs, and counted toward the genome size
                the location namespace has to cover.  See CONTIG_PADDING.

        Raises:
            SnapError: if the build fails, or runs out of location namespace
                at every size tried.
        """
        fasta = Path(fasta)
        self.index_name = fasta.stem
        index_dir = Path(self.index_name)
        self._lock_index(index_dir)

        # Require the marker, not just the directory: a build killed partway
        # leaves the directory behind, and trusting it yields a corrupt index.
        if (index_dir / _INDEX_MARKER).exists():
            logger.info("reusing snap index %s", index_dir.resolve())
            self.index_built = True
            return self.index_name

        # Logged at INFO, not DEBUG: an index build is the longest step in a
        # run and rebuilding one that should have been reused is invisible
        # otherwise -- the reuse branch above returns without a word, so a log
        # that jumps straight to the -locationSize sweep looks identical
        # whether the index was missing or was deleted between runs.
        logger.info(
            "no snap index at %s (looked for %s); building one",
            index_dir.resolve(),
            _INDEX_MARKER,
        )

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
                f"-p{padding}",         # see CONTIG_PADDING
                "-locationSize", size,
            ]
            result = run(args, log_path=LOG_PATH)
            if result.ok:
                # SILENT_OPTION_REJECTION applies to indexing too: snap can
                # exit 0 having written nothing, and mapping against an index
                # that is not there dies with no more than snap's banner --
                # which is a much harder thing to read than this message.
                if not (index_dir / _INDEX_MARKER).exists():
                    raise SnapError(
                        f"snap-aligner exited successfully but wrote no index "
                        f"at {index_dir.resolve()} (no {_INDEX_MARKER}).\n"
                        f"command: {' '.join(str(a) for a in args)}\n"
                        f"{_tail(result.output)}"
                    )
                logger.info(
                    "snap index built at -locationSize %d: %s",
                    size,
                    index_dir.resolve(),
                )
                self.index_built = True
                return self.index_name

            last_error = result.stderr or result.stdout
            if any(p.search(last_error) for p in _LOCATION_SIZE_PATTERNS):
                # Clear the partial index before retrying. rmdir() will not
                # do -- snap leaves Genome, GenomeIndexHash and OverflowTable
                # behind, and rmdir only removes empty dirs.
                #
                # Never a complete one, though: this deletes the longest
                # piece of work in a run, and an index carrying the marker
                # belongs to whoever built it.  The check at the top of
                # build_index should already have returned in that case, so
                # reaching here means something is wrong -- say so rather
                # than destroying the index and rebuilding it next run.
                if (index_dir / _INDEX_MARKER).exists():
                    raise SnapError(
                        f"snap says the location size is too small, but "
                        f"{index_dir.resolve()} holds a complete index. "
                        f"Refusing to delete it. Remove it by hand if it is "
                        f"stale, or pass --location-size to skip the sweep."
                        f"\n{_tail(last_error)}"
                    )
                if index_dir.is_dir():
                    shutil.rmtree(index_dir, ignore_errors=True)
                if size == sizes[0]:
                    # Once per build, not once per attempt: the figures do
                    # not change as the sweep climbs, and reading the FASTA
                    # again for each would not be free.
                    report = _padding_report(fasta, padding)
                    if report:
                        logger.warning("%s", report)
                if size == sizes[-1]:
                    break
                logger.warning(
                    "snap ran out of genome locations at -locationSize %d, "
                    "retrying at %d (a full rebuild, and a larger index for "
                    "the rest of the run)",
                    size,
                    size + 1,
                )
                continue
            raise SnapError(
                f"Failed to build snap index: "
                f"{_how_it_died(result.returncode)}. "
                f"Full output in {LOG_PATH.resolve()}\n{_tail(last_error)}"
            )

        hint = (
            " Every location size from "
            f"{sizes[0]} to {sizes[-1]} was too small; try a larger --seed-size."
            if len(sizes) > 1
            else f" Retry without --location-size to sweep {_LOCATION_SIZES.start}"
                 f"-{_LOCATION_SIZES.stop - 1}, or use a larger --seed-size."
        )
        raise SnapError(f"Failed to build snap index.{hint}\n{last_error}")

    # -- mapping ----------------------------------------------------------

    def build_paired_command(
        self,
        left,
        right,
        threads: int,
        output: str,
        max_seed_hits: int = DEFAULT_MAX_SEED_HITS,
        edit_distance: int = DEFAULT_EDIT_DISTANCE,
        extra_search_depth: int = DEFAULT_EXTRA_SEARCH_DEPTH,
        multi_edit_distance: int = DEFAULT_MULTI_EDIT_DISTANCE,
        max_alignments_per_pair: int = DEFAULT_MAX_ALIGNMENTS_PER_PAIR,
        max_alignments_per_contig: int | None = DEFAULT_MAX_ALIGNMENTS_PER_CONTIG,
        max_candidate_pool: int | None = None,
    ):
        """Assemble the ``snap-aligner paired`` command.

        See MULTI_ALIGNMENT_SETTINGS for why the -D/-om/-omax/-mpc defaults
        differ from the Ruby's, and MAX_CANDIDATE_POOL for why -mcp is gone.
        """
        args = [self.binary, "paired", self.index_name]
        for l, r in zip(str(left).split(","), str(right).split(",")):
            args += [l, r]
        args += [
            "-o", output,
            "-s", 0, 1000,        # min/max spacing between paired-read starts
            "-H", max_seed_hits,  # max seed hits in paired mode
            "-h", 2000,           # max seed hits when reverting to single
            "-d", edit_distance,  # max edit distance
            "-t", threads,
            "-b",                 # bind threads to cores
            "-M",                 # M-style CIGAR (now the default, kept explicit)
            "-D", extra_search_depth,       # must track -om
            "-om", multi_edit_distance,     # report multiple alignments
            "-omax", max_alignments_per_pair,
        ]
        if max_alignments_per_contig is not None:
            # Applied before -omax: best placement per candidate contig,
            # which is what the fragment assignment actually consumes.
            args += ["-mpc", max_alignments_per_contig]
        if max_candidate_pool is not None:
            args += ["-mcp", max_candidate_pool]
        return args

    def map_reads(
        self,
        left,
        right,
        threads: int = 8,
        output=None,
        max_seed_hits: int = DEFAULT_MAX_SEED_HITS,
        edit_distance: int = DEFAULT_EDIT_DISTANCE,
        extra_search_depth: int = DEFAULT_EXTRA_SEARCH_DEPTH,
        multi_edit_distance: int = DEFAULT_MULTI_EDIT_DISTANCE,
        max_alignments_per_pair: int = DEFAULT_MAX_ALIGNMENTS_PER_PAIR,
        max_alignments_per_contig: int | None = DEFAULT_MAX_ALIGNMENTS_PER_CONTIG,
        max_candidate_pool: int | None = None,
    ) -> str:
        """Map paired reads, returning the path to the BAM.

        The BAM is left in **read order**, not coordinate sorted: both
        :func:`~pytransrate.assign.assign_fragments` and
        :func:`~pytransrate.bam_metrics.estimate_realistic_distance` require
        mates and multi-mappings to be adjacent.
        """
        if not self.index_built:
            raise SnapError("Index not built")

        lbase = Path(str(left).split(",")[0]).name
        rbase = Path(str(right).split(",")[0]).name
        index = Path(self.index_name).name
        self.bam = str(Path(output or f"{lbase}.{rbase}.{index}.bam").resolve())
        self._read_count_file = f"{lbase}-{rbase}-read_count.txt"

        # Logged, both ways, for the same reason the index build is (see
        # build_index): reuse used to be silent, so a log that jumps straight
        # into snap's output looked identical whether the BAM was missing,
        # was unusable, or was reused -- and this is the step that costs
        # hours on a real library.
        if os.path.exists(self.bam):
            size_gb = os.path.getsize(self.bam) / 1e9
            if alignment_is_done(self.bam):
                logger.info(
                    "reusing existing BAM (%.1f GB), %s says mapping finished: %s",
                    size_gb,
                    align_done_path(self.bam).name,
                    self.bam,
                )
                self._load_read_count(left)
                return self.bam
            if bam_is_complete(self.bam):
                # Predates ALIGN_DONE, or the marker was removed by hand.
                # The EOF marker is the same evidence the check has always
                # used, so this stays a reuse rather than hours re-spent.
                logger.info(
                    "reusing existing BAM (%.1f GB): it has no %s but ends "
                    "with the BGZF end-of-file marker, so mapping finished: %s",
                    size_gb,
                    ALIGN_DONE_SUFFIX,
                    self.bam,
                )
                self._write_align_done(self.bam, note="marker added on reuse")
                self._load_read_count(left)
                return self.bam
            # See ALIGN_DONE: unusable for metrics, but not ours to destroy.
            kept = Path(str(self.bam) + PARTIAL_SUFFIX)
            try:
                os.replace(self.bam, kept)
                logger.warning(
                    "%s (%.1f GB) stops before the BGZF end-of-file marker, "
                    "so mapping did not finish -- a killed run leaves exactly "
                    "this. It is NOT being overwritten: moved to %s. Delete "
                    "it when you no longer need it (rm %s). Mapping again.",
                    self.bam, size_gb, kept.name, kept,
                )
            except OSError as exc:
                logger.warning(
                    "%s (%.1f GB) is a partial BAM from a run that did not "
                    "finish, and it could not be moved aside (%s), so snap "
                    "will overwrite it. Mapping again.",
                    self.bam, size_gb, exc,
                )
        else:
            logger.info("no BAM at %s; mapping", self.bam)

        args = self.build_paired_command(
            left, right, threads, self.bam,
            max_seed_hits=max_seed_hits,
            edit_distance=edit_distance,
            extra_search_depth=extra_search_depth,
            multi_edit_distance=multi_edit_distance,
            max_alignments_per_pair=max_alignments_per_pair,
            max_alignments_per_contig=max_alignments_per_contig,
            max_candidate_pool=max_candidate_pool,
        )
        result = run(args, log_path=LOG_PATH)
        self._save_read_count(result.output)

        # SILENT_OPTION_REJECTION: check this before the exit code, because
        # snap returns 0 in this case.
        bad_options = _BAD_OPTIONS.search(result.output)
        if bad_options:
            raise SnapError(
                "snap-aligner rejected an option and produced no alignments: "
                f"{bad_options.group(1).strip()}\n"
                f"command: {' '.join(str(a) for a in args)}"
            )

        if not result.ok:
            match = _UNMATCHED_IDS.search(result.output)
            if match:
                raise SnapError(
                    "snap found unmatched read IDs in the input fastq files.\n"
                    f"Left files contained read id\n{match.group(1).strip()}\n"
                    f"and right files contained read id\n{match.group(2).strip()}\n"
                    "at the same position in the file."
                )
            raise SnapError(
                f"snap failed: {_how_it_died(result.returncode)}. "
                f"Full output in {LOG_PATH.resolve()}\n{_tail(result.output)}"
            )

        # snap can exit 0 having written nothing; see SILENT_OPTION_REJECTION.
        if not os.path.exists(self.bam) or os.path.getsize(self.bam) == 0:
            raise SnapError(
                f"snap-aligner exited successfully but produced no alignments "
                f"at {self.bam}\ncommand: "
                f"{' '.join(str(a) for a in args)}\n{_tail(result.output)}"
            )

        # Last thing, deliberately: everything above is a way for a snap that
        # returned 0 to still have produced nothing usable, and the marker
        # must mean all of them passed. See ALIGN_DONE.
        self._write_align_done(self.bam, command=args)
        return self.bam

    def _write_align_done(self, bam, command=None, note=None) -> None:
        """Record that mapping finished. See ALIGN_DONE.

        Never fatal: the marker is an optimisation for the next run, and a
        run that mapped successfully must not be failed by a read-only
        output directory at the very end of it.
        """
        try:
            lines = [
                f"finished: {datetime.datetime.now().isoformat(timespec='seconds')}",
                f"bam: {bam}",
                f"bytes: {os.path.getsize(bam)}",
                f"fragments: {self.read_count}",
            ]
            if command is not None:
                lines.append("command: " + " ".join(str(a) for a in command))
            if note is not None:
                lines.append(f"note: {note}")
            align_done_path(bam).write_text("\n".join(lines) + "\n")
        except OSError as exc:
            logger.warning(
                "mapping finished but %s could not be written (%s); the next "
                "run falls back to checking the BAM's end-of-file marker",
                align_done_path(bam), exc,
            )

    # -- read counting ----------------------------------------------------

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

        logger.info(
            "no saved read count at %s; counting the reads in %s instead, "
            "which reads the whole library",
            self._read_count_file,
            reads,
        )
        for path in str(reads).split(","):
            try:
                # open_binary, not open: counting the lines of a gzip stream
                # counts compressed data and gives a meaningless figure.
                with open_binary(path) as handle:
                    lines = sum(1 for _ in handle)
                self.read_count += lines // 4
            except (OSError, EOFError):  # EOFError: a truncated gzip file
                logger.warning("couldn't count reads in %s", path)

        if self._read_count_file:
            with open(self._read_count_file, "w") as handle:
                handle.write(f"{self.read_count}\n")
