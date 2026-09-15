"""Tests for the external-tool wrappers and the CLI.

Command construction is asserted directly rather than by running the tools,
so these pass without snap-aligner or salmon installed. The end-to-end run
against the real binaries lives in tests/test_pipeline.py.
"""

from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path

import pytest

from pytransrate import cmd
from pytransrate.cmd import CommandError, CommandResult, run, which
from pytransrate.mapper import Snap
from pytransrate.quantify import Salmon, SalmonError, load_expression
from pytransrate.read_metrics import get_read_length


# ---------------------------------------------------------------------------
# cmd
# ---------------------------------------------------------------------------


def test_run_captures_stdout():
    result = run(["echo", "hello"])
    assert result.ok
    assert result.stdout.strip() == "hello"


def test_run_coerces_arguments_to_strings():
    result = run(["echo", 42])
    assert result.stdout.strip() == "42"


def test_run_reports_failure_without_raising():
    result = run(["false"])
    assert not result.ok
    with pytest.raises(CommandError):
        result.check()


def test_run_does_not_use_a_shell():
    """Arguments are passed through, not interpreted."""
    result = run(["echo", "$HOME; rm -rf /"])
    assert "$HOME" in result.stdout


def test_run_writes_stdout_to_a_file(tmp_path):
    out = tmp_path / "o.txt"
    result = run(["echo", "written"], stdout_path=str(out))
    assert result.ok
    assert out.read_text().strip() == "written"


def test_which_raises_for_a_missing_binary():
    with pytest.raises(CommandError, match="could not find"):
        which("definitely-not-a-real-binary-xyz")


def test_command_result_str_is_the_command():
    assert str(CommandResult(args=["a", "b"], returncode=0)) == "a b"


# ---------------------------------------------------------------------------
# snap-aligner command construction
# ---------------------------------------------------------------------------


@pytest.fixture
def snap():
    obj = Snap.__new__(Snap)          # bypass the PATH lookup
    obj.binary = "snap-aligner"
    obj.index_name = "assembly"
    obj.index_built = True
    obj.bam = None
    obj.read_count = 0
    obj._read_count_file = None
    return obj


def test_paired_command_keeps_the_ruby_flag_set(snap):
    """The same flags as the Ruby; several values differ, see
    MULTI_ALIGNMENT_SETTINGS."""
    args = [str(a) for a in snap.build_paired_command("l.fq", "r.fq", 8, "o.bam")]
    joined = " ".join(args)
    for flag in ["-s", "-H", "-h", "-d", "-t", "-b", "-M", "-D", "-om",
                 "-omax", "-o"]:
        assert flag in args, flag
    assert "-s 0 1000" in joined
    assert args[:3] == ["snap-aligner", "paired", "assembly"]


def test_paired_command_interleaves_multiple_read_files(snap):
    args = [str(a) for a in
            snap.build_paired_command("a.fq,b.fq", "x.fq,y.fq", 4, "o.bam")]
    idx = args.index("a.fq")
    assert args[idx:idx + 4] == ["a.fq", "x.fq", "b.fq", "y.fq"]


def test_map_reads_requires_an_index(snap):
    snap.index_built = False
    with pytest.raises(Exception, match="Index not built"):
        snap.map_reads("l.fq", "r.fq")


def test_read_count_is_halved_from_snaps_summary(snap):
    snap._read_count_file = None
    snap._save_read_count(
        "Total Reads    Aligned    Unaligned   x   y   z\n"
        "1,000          900        100         1   2   3\n"
    )
    assert snap.read_count == 500


# ---------------------------------------------------------------------------
# salmon command construction and quant.sf parsing
# ---------------------------------------------------------------------------


@pytest.fixture
def salmon():
    obj = Salmon.__new__(Salmon)
    obj.binary = "salmon"
    return obj


def test_salmon_command_uses_modern_flags(salmon):
    args = [str(a) for a in salmon.build_command("a.fa", "in.bam", 8, "out")]
    joined = " ".join(args)
    assert "--alignments in.bam" in joined
    assert "--targets a.fa" in joined
    assert "--libType A" in joined
    assert "--errorModel" in args


def test_salmon_command_omits_removed_flags(salmon):
    """--useErrorModel is a hard error in salmon 2.x; --sampleOut is inert."""
    args = [str(a) for a in salmon.build_command("a.fa", "in.bam", 8, "out")]
    assert "--useErrorModel" not in args
    assert "--sampleOut" not in args
    assert "--sampleUnaligned" not in args


def test_salmon_error_model_can_be_disabled(salmon):
    args = [str(a) for a in
            salmon.build_command("a.fa", "in.bam", 8, "out", error_model=False)]
    assert "--errorModel" not in args


def _quant(tmp_path, rows):
    path = tmp_path / "quant.sf"
    path.write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        + "".join("\t".join(str(c) for c in r) + "\n" for r in rows)
    )
    return path


def test_load_expression_parses_quant_sf(tmp_path):
    path = _quant(tmp_path, [["tx0", 1000, 850.5, 123.4, 42.7]])
    expression = load_expression(path)
    assert expression["tx0"]["eff_len"] == 850     # truncated, as in the Ruby
    assert expression["tx0"]["tpm"] == pytest.approx(123.4)
    assert expression["tx0"]["eff_count"] == pytest.approx(42.7)


def test_load_expression_rejects_the_wrong_column_count(tmp_path):
    path = tmp_path / "quant.sf"
    path.write_text("Name\tLength\tTPM\n" + "tx0\t100\t1.0\n")
    with pytest.raises(SalmonError, match="5 columns"):
        load_expression(path)


def test_load_expression_on_an_empty_table(tmp_path):
    assert load_expression(_quant(tmp_path, [])) == {}


# ---------------------------------------------------------------------------
# read length
# ---------------------------------------------------------------------------


def test_get_read_length_takes_the_maximum(tmp_path):
    fq = tmp_path / "r.fq"
    fq.write_text(
        "@a\n" + "A" * 75 + "\n+\n" + "I" * 75 + "\n"
        "@b\n" + "A" * 100 + "\n+\n" + "I" * 100 + "\n"
    )
    assert get_read_length(str(fq)) == 100


def test_get_read_length_reads_gzip(tmp_path):
    fq = tmp_path / "r.fq.gz"
    with gzip.open(fq, "wt") as handle:
        handle.write("@a\n" + "A" * 60 + "\n+\n" + "I" * 60 + "\n")
    assert get_read_length(str(fq)) == 60


def test_get_read_length_uses_the_first_file_only(tmp_path):
    a = tmp_path / "a.fq"
    a.write_text("@a\n" + "A" * 50 + "\n+\n" + "I" * 50 + "\n")
    b = tmp_path / "b.fq"
    b.write_text("@b\n" + "A" * 200 + "\n+\n" + "I" * 200 + "\n")
    assert get_read_length(f"{a},{b}") == 50


# ---------------------------------------------------------------------------
# CLI argument validation
# ---------------------------------------------------------------------------


def _args(**kwargs):
    from pytransrate.cli import build_parser

    argv = []
    for key, value in kwargs.items():
        flag = "--" + key.replace("_", "-")
        if value is True:
            argv.append(flag)
        elif value is not None:
            argv += [flag, str(value)]
    return build_parser().parse_args(argv)


def test_reference_is_rejected_with_an_explanation(tmp_path):
    from pytransrate.cli import check_arguments

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    args = _args(assembly=str(fasta), reference=str(fasta))
    with pytest.raises(CommandError, match="not implemented"):
        check_arguments(args)


def test_missing_assembly_is_rejected():
    from pytransrate.cli import check_arguments

    with pytest.raises(CommandError, match="does not exist"):
        check_arguments(_args(assembly="/nope/missing.fa"))


def test_duplicate_assemblies_are_rejected(tmp_path):
    from pytransrate.cli import check_arguments

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    args = _args(assembly=f"{fasta},{fasta}")
    with pytest.raises(CommandError, match="more than once"):
        check_arguments(args)


def test_left_without_right_is_rejected(tmp_path):
    from pytransrate.cli import check_arguments

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    reads = tmp_path / "r.fq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    args = _args(assembly=str(fasta), left=str(reads))
    with pytest.raises(CommandError, match="together"):
        check_arguments(args)


def test_mismatched_read_file_counts_are_rejected(tmp_path):
    from pytransrate.cli import check_arguments

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    r1 = tmp_path / "r1.fq"
    r1.write_text("@a\nACGT\n+\nIIII\n")
    args = _args(assembly=str(fasta), left=f"{r1},{r1}", right=str(r1))
    with pytest.raises(CommandError, match="same number"):
        check_arguments(args)


def test_read_paths_are_absolutised(tmp_path, monkeypatch):
    """They must resolve before the run chdirs into the output directory."""
    from pytransrate.cli import check_arguments

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    reads = tmp_path / "r.fq"
    reads.write_text("@a\nACGT\n+\nIIII\n")

    monkeypatch.chdir(tmp_path)
    args = _args(assembly="a.fa", left="r.fq", right="r.fq")
    check_arguments(args)
    assert args.left.startswith("/")
    assert args.right.startswith("/")


# ---------------------------------------------------------------------------
# snap index: locationSize sweep and cleanup
# ---------------------------------------------------------------------------


class _Result:
    def __init__(self, ok, stderr=""):
        self.ok = ok
        self.returncode = 0 if ok else 1
        self.stderr = stderr
        self.stdout = ""

    @property
    def output(self):
        return self.stdout + self.stderr


_OVERFLOW = (
    "Ran out of overflow table namespace. This genome cannot be indexed "
    "with this seed and location size.  Increase at least one.\n"
)

# The other three ways snap 2.0.5 says the location size is too small, copied
# verbatim from GenomeIndex.cpp (double spaces included). Each must drive the
# sweep: matching only some of them leaves the user to set --location-size by
# hand, which is what pytransrate is meant to spare them.
_TOO_BIG = (
    "Welcome to SNAP version 2.0.5.\n\n"
    "Genome is too big for 4 byte genome locations.  Specify a larger "
    "location size with -locationSize\n"
    "SNAP exited with exit code 1 from line 569 of file "
    "SNAPLib/GenomeIndex.cpp\n"
)
_TOO_MANY_ENTRIES = (
    "Trying to use too many overflow entries.  To index this genome, you "
    "either need a larger seed size or a larger location size.\n"
)
_NO_ADDRESS_SPACE = (
    "Not enough address space to index this genome with this seed size.  "
    "Try a larger seed or location size.\n"
)


def _fake_run(monkeypatch, outcomes, calls):
    import pytransrate.mapper as mapper

    def fake(args, **kwargs):
        args = [str(a) for a in args]
        calls.append(args)
        # Materialise what snap leaves behind, so cleanup is exercised.
        index_dir = Path(args[3])
        index_dir.mkdir(exist_ok=True)
        result = outcomes.pop(0)
        for name in ("Genome", "GenomeIndexHash", "OverflowTable"):
            (index_dir / name).write_text("x")
        if result.ok:
            (index_dir / "GenomeIndex").write_text("x")
        return result

    monkeypatch.setattr(mapper, "run", fake)


def _snap():
    obj = Snap.__new__(Snap)
    obj.binary = "snap-aligner"
    obj.index_name = None
    obj.index_built = False
    obj.bam = None
    obj.read_count = 0
    obj._read_count_file = None
    obj._lock = None
    return obj


def _location_sizes(calls):
    return [c[c.index("-locationSize") + 1] for c in calls]


def test_index_starts_at_location_size_four(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(True)], calls)
    _snap().build_index("a.fa")
    assert _location_sizes(calls) == ["4"]


def test_index_steps_up_on_overflow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(
        monkeypatch,
        [_Result(False, _OVERFLOW), _Result(False, _OVERFLOW), _Result(True)],
        calls,
    )
    _snap().build_index("a.fa")
    assert _location_sizes(calls) == ["4", "5", "6"]


@pytest.mark.parametrize(
    "message",
    [_TOO_BIG, _TOO_MANY_ENTRIES, _NO_ADDRESS_SPACE],
    ids=["too_big", "too_many_entries", "no_address_space"],
)
def test_every_location_size_message_drives_the_sweep(
    tmp_path, monkeypatch, message
):
    """Regression: only two of snap's four messages used to be recognised.

    A large merged assembly hits "Genome is too big for 4 byte genome
    locations" first, and that one fell through to a hard failure -- the
    sweep never ran and the build died at -locationSize 4.
    """
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(False, message), _Result(True)], calls)
    _snap().build_index("a.fa")
    assert _location_sizes(calls) == ["4", "5"]


def test_sweep_message_on_stdout_is_recognised(tmp_path, monkeypatch):
    """snap writes these to stderr, but run() falls back to stdout."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    failed = _Result(False)
    failed.stdout = _TOO_BIG
    _fake_run(monkeypatch, [failed, _Result(True)], calls)
    _snap().build_index("a.fa")
    assert _location_sizes(calls) == ["4", "5"]


def test_partial_index_is_removed_between_attempts(tmp_path, monkeypatch):
    """rmdir() could not do this -- snap leaves four files behind."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    seen = []

    import pytransrate.mapper as mapper

    outcomes = [_Result(False, _OVERFLOW), _Result(True)]

    def fake(args, **kwargs):
        args = [str(a) for a in args]
        index_dir = Path(args[3])
        seen.append(sorted(p.name for p in index_dir.iterdir())
                    if index_dir.is_dir() else [])
        index_dir.mkdir(exist_ok=True)
        for name in ("Genome", "GenomeIndexHash", "OverflowTable"):
            (index_dir / name).write_text("x")
        result = outcomes.pop(0)
        if result.ok:
            (index_dir / "GenomeIndex").write_text("x")
        return result

    monkeypatch.setattr(mapper, "run", fake)
    _snap().build_index("a.fa")
    # The second attempt must start from a clean directory.
    assert seen[1] == []


def test_index_raises_when_every_size_overflows(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(False, _OVERFLOW)] * 5, calls)
    with pytest.raises(Exception, match="seed-size"):
        _snap().build_index("a.fa")
    assert _location_sizes(calls) == ["4", "5", "6", "7", "8"]


def test_explicit_location_size_skips_the_sweep(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(True)], calls)
    _snap().build_index("a.fa", location_size=6)
    assert _location_sizes(calls) == ["6"]


def test_explicit_location_size_does_not_retry(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(False, _OVERFLOW)], calls)
    with pytest.raises(Exception, match="location-size"):
        _snap().build_index("a.fa", location_size=6)
    assert _location_sizes(calls) == ["6"]


def test_location_size_out_of_range_is_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    with pytest.raises(Exception, match="between 4 and 8"):
        _snap().build_index("a.fa", location_size=9)


def test_non_overflow_failure_does_not_sweep(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(False, "disk on fire")], calls)
    with pytest.raises(Exception, match="disk on fire"):
        _snap().build_index("a.fa")
    assert len(calls) == 1


def test_complete_index_is_reused(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    index = tmp_path / "a"
    index.mkdir()
    (index / "GenomeIndex").write_text("x")
    calls = []
    _fake_run(monkeypatch, [], calls)
    _snap().build_index("a.fa")
    assert calls == []


def test_partial_index_is_not_mistaken_for_a_complete_one(tmp_path, monkeypatch):
    """A directory without the marker must be rebuilt, not trusted."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    index = tmp_path / "a"
    index.mkdir()
    (index / "Genome").write_text("x")     # partial: no GenomeIndex marker
    calls = []
    _fake_run(monkeypatch, [_Result(True)], calls)
    _snap().build_index("a.fa")
    assert len(calls) == 1


def _padding(calls):
    return [a for c in calls for a in c if a.startswith("-p")]


def test_index_pads_below_snaps_own_default(tmp_path, monkeypatch):
    """snap's 2000 is sized for genomes; see CONTIG_PADDING in mapper."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(True)], calls)
    _snap().build_index("a.fa")
    assert _padding(calls) == ["-p1000"]


def test_padding_is_attached_to_the_flag_without_a_space(tmp_path, monkeypatch):
    """snap parses -p with atoi(argv[n] + 2), so a separate argument is lost."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(True)], calls)
    _snap().build_index("a.fa", padding=250)
    assert _padding(calls) == ["-p250"]
    assert "250" not in calls[0][calls[0].index("-p250") + 1:]


def test_padding_survives_the_location_size_sweep(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(
        monkeypatch,
        [_Result(False, "Genome is too big for 4 byte genome locations"),
         _Result(True)],
        calls,
    )
    _snap().build_index("a.fa")
    assert _location_sizes(calls) == ["4", "5"]
    assert _padding(calls) == ["-p1000", "-p1000"]


def test_cli_exposes_both_index_knobs():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(
        ["-a", "x.fa", "--location-size", "6", "--seed-size", "25"]
    )
    assert args.location_size == 6
    assert args.seed_size == 25


def test_cli_exposes_padding():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(["-a", "x.fa", "--padding", "1500"])
    assert args.padding == 1500


def test_cli_defaults_leave_the_sweep_enabled():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(["-a", "x.fa"])
    assert args.location_size is None
    assert args.seed_size == 23
    assert args.padding == 1000


# ---------------------------------------------------------------------------
# snap paired: tunables
# ---------------------------------------------------------------------------


def test_mcp_is_not_passed_by_default(snap):
    """MAX_CANDIDATE_POOL: the Ruby's value overflowed snap's atoi()."""
    args = [str(a) for a in snap.build_paired_command("l.fq", "r.fq", 8, "o.bam")]
    assert "-mcp" not in args


def test_mcp_is_passed_when_asked_for(snap):
    args = [str(a) for a in snap.build_paired_command(
        "l.fq", "r.fq", 8, "o.bam", max_candidate_pool=1_000_000)]
    assert args[args.index("-mcp") + 1] == "1000000"


def test_seed_hits_and_edit_distance_are_tunable(snap):
    args = [str(a) for a in snap.build_paired_command(
        "l.fq", "r.fq", 8, "o.bam", max_seed_hits=4000, edit_distance=20)]
    assert args[args.index("-H") + 1] == "4000"
    assert args[args.index("-d") + 1] == "20"


def test_paired_defaults_are_the_validated_configuration(snap):
    """MULTI_ALIGNMENT_SETTINGS: the flags that run on real ORP data."""
    args = [str(a) for a in snap.build_paired_command("l.fq", "r.fq", 8, "o.bam")]
    assert args[args.index("-H") + 1] == "4000"
    assert args[args.index("-d") + 1] == "30"
    assert args[args.index("-D") + 1] == "2"
    assert args[args.index("-om") + 1] == "2"
    assert args[args.index("-omax") + 1] == "10"
    assert args[args.index("-mpc") + 1] == "1"


def test_extra_search_depth_is_at_least_the_multi_edit_distance(snap):
    """-om needs -D to search that far; the defaults must not violate it."""
    args = [str(a) for a in snap.build_paired_command("l.fq", "r.fq", 8, "o.bam")]
    assert int(args[args.index("-D") + 1]) >= int(args[args.index("-om") + 1])


def test_multi_alignment_flags_are_tunable(snap):
    args = [str(a) for a in snap.build_paired_command(
        "l.fq", "r.fq", 8, "o.bam", extra_search_depth=5,
        multi_edit_distance=5, max_alignments_per_pair=3)]
    assert args[args.index("-D") + 1] == "5"
    assert args[args.index("-om") + 1] == "5"
    assert args[args.index("-omax") + 1] == "3"


def test_mpc_can_be_disabled(snap):
    args = [str(a) for a in snap.build_paired_command(
        "l.fq", "r.fq", 8, "o.bam", max_alignments_per_contig=None)]
    assert "-mpc" not in args


def test_multi_alignment_is_never_silently_disabled(snap):
    """-om is what produces the alternates the assignment step needs."""
    args = [str(a) for a in snap.build_paired_command("l.fq", "r.fq", 8, "o.bam")]
    assert "-om" in args
    assert int(args[args.index("-om") + 1]) >= 1


# ---------------------------------------------------------------------------
# SILENT_OPTION_REJECTION
# ---------------------------------------------------------------------------


def _mapping_snap(tmp_path, monkeypatch, result, write_bam=True, size=10):
    import pytransrate.mapper as mapper

    def fake(args, **kwargs):
        if write_bam:
            out = args[args.index("-o") + 1]
            Path(out).write_bytes(b"x" * size)
        return result

    monkeypatch.setattr(mapper, "run", fake)
    obj = Snap.__new__(Snap)
    obj.binary = "snap-aligner"
    obj.index_name = "assembly"
    obj.index_built = True
    obj.bam = None
    obj.read_count = 0
    obj._read_count_file = None
    obj._lock = None
    monkeypatch.chdir(tmp_path)
    return obj


def test_rejected_option_raises_despite_exit_zero(tmp_path, monkeypatch):
    """snap prints the complaint, writes nothing, and exits 0."""
    result = _Result(True)
    result.stdout = "Didn't understand options starting at -om 2 -omax 10\nUsage:"
    snap = _mapping_snap(tmp_path, monkeypatch, result, write_bam=False)
    with pytest.raises(Exception, match="rejected an option"):
        snap.map_reads("l.fq", "r.fq")


def test_missing_bam_raises_despite_exit_zero(tmp_path, monkeypatch):
    snap = _mapping_snap(tmp_path, monkeypatch, _Result(True), write_bam=False)
    with pytest.raises(Exception, match="produced no alignments"):
        snap.map_reads("l.fq", "r.fq")


def test_empty_bam_raises_despite_exit_zero(tmp_path, monkeypatch):
    snap = _mapping_snap(tmp_path, monkeypatch, _Result(True), size=0)
    with pytest.raises(Exception, match="produced no alignments"):
        snap.map_reads("l.fq", "r.fq")


def test_successful_mapping_returns_the_bam(tmp_path, monkeypatch):
    snap = _mapping_snap(tmp_path, monkeypatch, _Result(True))
    bam = snap.map_reads("l.fq", "r.fq")
    assert Path(bam).exists()


def test_cli_exposes_the_multi_alignment_flags():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(
        ["-a", "x.fa", "--extra-search-depth", "5", "--multi-edit-distance", "5",
         "--max-alignments-per-pair", "3", "--max-alignments-per-contig", "2"]
    )
    assert args.extra_search_depth == 5
    assert args.multi_edit_distance == 5
    assert args.max_alignments_per_pair == 3
    assert args.max_alignments_per_contig == 2


def test_cli_multi_alignment_defaults():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(["-a", "x.fa"])
    assert (args.max_seed_hits, args.extra_search_depth,
            args.multi_edit_distance, args.max_alignments_per_pair,
            args.max_alignments_per_contig) == (4000, 2, 2, 10, 1)


def test_cli_exposes_the_paired_tunables():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(
        ["-a", "x.fa", "--max-seed-hits", "4000", "--edit-distance", "20",
         "--max-candidate-pool", "1000000"]
    )
    assert args.max_seed_hits == 4000
    assert args.edit_distance == 20
    assert args.max_candidate_pool == 1_000_000


def test_cli_default_seed_hits_is_snaps_own():
    from pytransrate.cli import build_parser

    assert build_parser().parse_args(["-a", "x.fa"]).max_seed_hits == 4000


def test_cli_defaults_omit_mcp():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(["-a", "x.fa"])
    assert args.max_candidate_pool is None
    assert args.edit_distance == 30


# ---------------------------------------------------------------------------
# banner
# ---------------------------------------------------------------------------


def test_banner_always_carries_the_version():
    from pytransrate import __version__
    from pytransrate.banner import banner

    for width in (40, 60, 80, 200):
        assert __version__ in banner(width=width), width


def test_banner_names_the_tool_and_what_it_does():
    from pytransrate.banner import TAGLINE, banner

    wide = banner(width=100)
    assert TAGLINE in wide
    assert "transcriptome" in TAGLINE

    narrow = banner(width=50)
    assert "pytransrate" in narrow
    assert TAGLINE in narrow


def test_banner_falls_back_when_the_terminal_is_narrow():
    from pytransrate.banner import banner

    narrow = banner(width=50)
    assert max(len(line) for line in narrow.split("\n")) <= 60
    # The wordmark would wrap and look broken at this width.
    assert "░" not in narrow


def test_banner_colour_is_opt_in():
    from pytransrate.banner import banner

    assert "\033[" not in banner(colour=False, width=100)
    assert "\033[" in banner(colour=True, width=100)


def test_banner_respects_no_color(monkeypatch):
    import io

    from pytransrate.banner import _use_colour

    tty = io.StringIO()
    tty.isatty = lambda: True
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm")
    assert _use_colour(tty)
    monkeypatch.setenv("NO_COLOR", "1")
    assert not _use_colour(tty)


def test_banner_is_plain_when_not_a_terminal():
    import io

    from pytransrate.banner import _use_colour

    assert not _use_colour(io.StringIO())


def test_banner_goes_to_stderr_leaving_stdout_clean(capsys):
    from pytransrate.banner import print_banner

    print_banner()
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "pytransrate" in captured.err or "░" in captured.err


def test_no_banner_flag_exists():
    from pytransrate.cli import build_parser

    assert build_parser().parse_args(["-a", "x.fa"]).no_banner is False
    assert build_parser().parse_args(["-a", "x.fa", "--no-banner"]).no_banner


# ---------------------------------------------------------------------------
# REPORTED_METRICS
# ---------------------------------------------------------------------------


def test_the_run_report_goes_to_stdout(capsys):
    """The banner is decoration and stays on stderr; the numbers are data."""
    import logging

    from pytransrate.cli import configure_logging

    root = logging.getLogger()
    saved = root.handlers[:]
    root.handlers = []
    try:
        configure_logging("info")
        logging.getLogger("pytransrate").info("a metric")
    finally:
        for handler in root.handlers:
            handler.close()
        root.handlers = saved

    captured = capsys.readouterr()
    assert "a metric" in captured.out
    assert "a metric" not in captured.err


def test_the_command_is_recorded_verbatim():
    """A log that does not say which settings produced it cannot be matched
    to a run, which is the whole point of comparing runs."""
    from pytransrate.cli import invocation

    line = invocation(["-a", "asm.fa", "-o", "out dir", "--max-alignments-per-contig", "0"])

    assert "--max-alignments-per-contig 0" in line
    assert "'out dir'" in line  # quoted, so it pastes back


def test_contig_report_is_min_max_and_n50(caplog):
    """Three lines on purpose: the rest of the length distribution is in the
    CSV and is not what anyone watches a run for."""
    import logging

    from pytransrate.cli import REPORTED_CONTIG_KEYS, log_metrics

    stats = {
        "n_seqs": 10, "smallest": 201, "largest": 12000, "n_bases": 50000,
        "mean_len": 5000.0, "n50": 9000, "n90": 300, "gc": 0.42,
    }
    with caplog.at_level(logging.INFO, logger="pytransrate"):
        log_metrics("contig metrics", stats, REPORTED_CONTIG_KEYS)

    reported = [record.getMessage() for record in caplog.records]
    assert any("min contig length" in line and "201" in line for line in reported)
    assert any("max contig length" in line and "12,000" in line for line in reported)
    assert any("N50" in line and "9,000" in line for line in reported)
    assert not any("n90" in line or "mean_len" in line for line in reported)


def test_mapping_report_covers_every_read_stat(caplog):
    """Complete on purpose: these are what move when the aligner settings
    move, so omitting any would force a trip to the CSV."""
    import logging

    from pytransrate.cli import log_metrics
    from pytransrate.output import READ_STATS_KEYS

    stats = {key: 1 for key in READ_STATS_KEYS}
    with caplog.at_level(logging.INFO, logger="pytransrate"):
        log_metrics("mapping metrics", stats, READ_STATS_KEYS)

    reported = "\n".join(record.getMessage() for record in caplog.records)
    for key in READ_STATS_KEYS:
        assert key in reported


def test_reported_values_match_the_csv_rounding():
    """A number read off the log must be the number in assemblies.csv,
    which write_assemblies_csv rounds to 5 places."""
    from pytransrate.cli import format_metric

    assert format_metric(0.8624231) == "0.86242"
    assert format_metric(28976658) == "28,976,658"
    assert format_metric(0) == "0"


# ---------------------------------------------------------------------------
# LIVE_LOGGING and INDEX_LOCK
# ---------------------------------------------------------------------------


def test_log_path_records_the_command_and_both_streams(tmp_path):
    log = tmp_path / "logs" / "snap.log"
    result = run(
        ["sh", "-c", "echo to-stdout; echo to-stderr >&2; exit 3"], log_path=log
    )

    assert result.returncode == 3
    assert "to-stdout" in result.output and "to-stderr" in result.output
    written = log.read_text()
    assert "to-stdout" in written and "to-stderr" in written
    assert written.startswith("$ sh -c ")


def test_log_survives_a_command_that_is_killed(tmp_path):
    """The whole point: an OOM-killed snap still explains itself.

    capture_output loses everything here, which is why the child is handed a
    file descriptor instead.
    """
    log = tmp_path / "logs" / "snap.log"
    result = run(["sh", "-c", "echo got-this-far; kill -9 $$"], log_path=log)

    assert not result.ok
    assert "got-this-far" in log.read_text()
    assert "got-this-far" in result.output


@pytest.fixture
def spawned_command(monkeypatch):
    """Record the argv _run_logged actually spawns, without spawning it."""
    recorded = []

    def fake_run(args, **kwargs):
        recorded.append(args)
        return subprocess.CompletedProcess(args, 0)

    monkeypatch.setattr(cmd.subprocess, "run", fake_run)
    cmd._line_buffered.cache_clear()
    yield recorded
    cmd._line_buffered.cache_clear()


def test_a_logged_child_is_line_buffered_so_a_crash_keeps_its_last_words(
    tmp_path, monkeypatch, spawned_command
):
    """LIVE_LOGGING: the file descriptor is useless if the child buffers.

    snap's account of a SIGFPE sat in its 4 KB stdio buffer and died with it,
    leaving a log holding only the banner snap had written to stderr.
    """
    monkeypatch.setattr(cmd.shutil, "which", lambda name: f"/usr/bin/{name}")

    result = run(["snap-aligner", "paired"], log_path=tmp_path / "snap.log")

    assert spawned_command[0] == [
        "/usr/bin/stdbuf", "-oL", "-eL", "snap-aligner", "paired",
    ]
    # ...but the command handed back to the user is one they can rerun.
    assert result.args == ["snap-aligner", "paired"]


def test_line_buffering_is_skipped_where_stdbuf_does_not_exist(
    tmp_path, monkeypatch, spawned_command
):
    """macOS ships no stdbuf; the run proceeds, merely as mute as before."""
    monkeypatch.setattr(cmd.shutil, "which", lambda name: None)

    run(["snap-aligner", "paired"], log_path=tmp_path / "snap.log")

    assert spawned_command[0] == ["snap-aligner", "paired"]


def test_log_appends_rather_than_replacing_earlier_runs(tmp_path):
    """build_index can run snap five times; each must stay readable."""
    log = tmp_path / "logs" / "snap.log"
    run(["echo", "first"], log_path=log)
    second = run(["echo", "second"], log_path=log)

    assert "first" in log.read_text()
    # ...but only this run's output comes back, or the read count would be
    # parsed out of a previous run's summary table.
    assert "first" not in second.output
    assert "second" in second.output


def test_a_second_run_cannot_delete_an_index_in_use(tmp_path, monkeypatch):
    """INDEX_LOCK: the loser is told to move, not left to trample."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    calls = []
    _fake_run(monkeypatch, [_Result(True), _Result(True)], calls)

    first = _snap()
    first.build_index("a.fa")

    with pytest.raises(Exception, match=f"in use by pid {os.getpid()}"):
        _snap().build_index("a.fa")

    # Whatever the second run wanted, it did not get as far as snap.
    assert len(calls) == 1


def test_closing_hands_the_index_to_the_next_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    _fake_run(monkeypatch, [_Result(True), _Result(True)], [])

    first = _snap()
    first.build_index("a.fa")
    first.close()

    second = _snap()
    assert second.build_index("a.fa") == "a"


def test_the_lock_sits_outside_the_index_directory(tmp_path, monkeypatch):
    """rmtree would otherwise delete the lock its own holder is holding,
    and the next process would lock a fresh inode unopposed."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    _fake_run(monkeypatch, [_Result(True)], [])

    _snap().build_index("a.fa")

    assert (tmp_path / "a.index.lock").exists()
    assert not (tmp_path / "a" / "a.index.lock").exists()


def test_exit_zero_without_an_index_is_not_success(tmp_path, monkeypatch):
    """snap can exit 0 having written nothing; mapping against the index
    that isn't there dies with nothing but snap's banner."""
    import pytransrate.mapper as mapper

    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    monkeypatch.setattr(mapper, "run", lambda args, **kw: _Result(True))

    with pytest.raises(Exception, match="wrote no index"):
        _snap().build_index("a.fa")


def test_a_complete_index_is_never_deleted_by_the_sweep(tmp_path, monkeypatch):
    """The sweep's cleanup must not take the longest piece of work in a run."""
    import pytransrate.mapper as mapper

    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")

    def fake(args, **kwargs):
        index_dir = Path(str(args[3]))
        index_dir.mkdir(exist_ok=True)
        # A complete index appears under the sweep's feet -- a concurrent
        # run, or a marker check that somehow did not fire.
        (index_dir / "GenomeIndex").write_text("x")
        return _Result(False, _OVERFLOW)

    monkeypatch.setattr(mapper, "run", fake)

    with pytest.raises(Exception, match="Refusing to delete"):
        _snap().build_index("a.fa")
    assert (tmp_path / "a" / "GenomeIndex").exists()


def test_reuse_and_rebuild_are_both_logged(tmp_path, monkeypatch, caplog):
    """A run that rebuilds an index it should have reused looks exactly like
    a first run unless the decision is logged."""
    import logging

    import pytransrate.mapper as mapper

    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c\nACGT\n")
    _fake_run(monkeypatch, [_Result(True)], [])

    with caplog.at_level(logging.INFO, logger="pytransrate"):
        builder = _snap()
        builder.build_index("a.fa")
        first = "\n".join(r.getMessage() for r in caplog.records)
        caplog.clear()
        builder.close()          # hand the lock over deterministically
        _snap().build_index("a.fa")
        again = "\n".join(r.getMessage() for r in caplog.records)

    assert "no snap index" in first and "built at -locationSize 4" in first
    assert "reusing snap index" in again


def test_location_size_failure_reports_the_padding_share(tmp_path, monkeypatch, caplog):
    """--padding, not --seed-size, is usually the term driving the ceiling."""
    import logging

    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.fa").write_text(">c1\nACGT\n>c2\nACGT\n>c3\nACGT\n")
    _fake_run(monkeypatch, [_Result(False, _TOO_BIG)] * 5, [])

    with caplog.at_level(logging.WARNING, logger="pytransrate"):
        with pytest.raises(Exception):
            _snap().build_index("a.fa")
    reported = "\n".join(r.getMessage() for r in caplog.records)

    assert "3 contigs x 1000 bp of padding" in reported
    assert "--padding" in reported
    # Once per build, however many sizes the sweep tries.
    assert reported.count("snap measures this genome") == 1


def test_contigs_are_counted_across_read_boundaries(tmp_path):
    """The count is chunked, so a '>' landing on a chunk edge must not be
    missed or doubled."""
    from pytransrate.mapper import _count_contigs

    fasta = tmp_path / "a.fa"
    fasta.write_text("".join(f">c{i}\n{'ACGT' * 500}\n" for i in range(2000)))
    assert _count_contigs(fasta) == 2000


def test_a_signal_death_is_named_not_left_as_a_negative_number(tmp_path, monkeypatch):
    """"-8" tells nobody anything; SIGFPE points straight at snap."""
    import pytransrate.mapper as mapper

    result = _Result(False)
    result.returncode = -8
    snap = _mapping_snap(tmp_path, monkeypatch, result)

    with pytest.raises(Exception, match="signal 8 .SIGFPE."):
        snap.map_reads("l.fq", "r.fq")


def test_a_kill_points_outward_rather_than_at_the_flags(tmp_path, monkeypatch):
    import pytransrate.mapper as mapper

    result = _Result(False)
    result.returncode = -9
    snap = _mapping_snap(tmp_path, monkeypatch, result)

    with pytest.raises(Exception, match="memory ceiling and wall clock"):
        snap.map_reads("l.fq", "r.fq")


def test_a_plain_failure_still_reads_as_an_exit_code(tmp_path, monkeypatch):
    result = _Result(False, "some snap complaint")
    result.returncode = 1
    snap = _mapping_snap(tmp_path, monkeypatch, result)

    with pytest.raises(Exception, match="exit 1"):
        snap.map_reads("l.fq", "r.fq")
