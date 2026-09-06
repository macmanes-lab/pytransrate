"""Tests for the external-tool wrappers and the CLI.

Command construction is asserted directly rather than by running the tools,
so these pass without snap-aligner or salmon installed. The end-to-end run
against the real binaries lives in tests/test_pipeline.py.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

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
        self.stderr = stderr
        self.stdout = ""


_OVERFLOW = (
    "Ran out of overflow table namespace. This genome cannot be indexed "
    "with this seed and location size.  Increase at least one.\n"
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


def test_cli_exposes_both_index_knobs():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(
        ["-a", "x.fa", "--location-size", "6", "--seed-size", "25"]
    )
    assert args.location_size == 6
    assert args.seed_size == 25


def test_cli_defaults_leave_the_sweep_enabled():
    from pytransrate.cli import build_parser

    args = build_parser().parse_args(["-a", "x.fa"])
    assert args.location_size is None
    assert args.seed_size == 23


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
