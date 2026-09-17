"""What survives a failure.

These cover the promise that a run which dies does not take the hours that
produced it along: the BAM is kept, a partial BAM is moved aside rather than
overwritten, and the marker that decides between the two is written only when
mapping really finished.  See ALIGN_DONE in pytransrate.mapper and
DEFERRED_BAM_DELETE in pytransrate.cli.

The aligner is faked throughout -- none of this needs snap, and the point is
the file lifecycle around it rather than anything snap does.
"""

import logging
from pathlib import Path

import pytest

from pytransrate.cmd import CommandError
from pytransrate.mapper import (
    Snap,
    align_done_path,
    alignment_is_done,
    bam_is_complete,
    _BGZF_EOF,
)


def _complete(path, size=200):
    """A BAM that ends the way a closed one does."""
    Path(path).write_bytes(b"BAM\1" + b"x" * size + _BGZF_EOF)


def _partial(path, size=200):
    """A BAM left by a run that was killed mid-write."""
    Path(path).write_bytes(b"BAM\1" + b"x" * size)


def _snap(tmp_path, monkeypatch, on_run=None):
    import pytransrate.mapper as mapper

    calls = []

    def fake(args, **kwargs):
        calls.append([str(a) for a in args])
        out = args[args.index("-o") + 1]
        if on_run is not None:
            on_run(out)
        else:
            _complete(out, size=10)

        class Result:
            ok = True
            returncode = 0
            stdout = "1,000 2 3 4 5 6\n"
            stderr = ""
            output = stdout

        return Result()

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
    return obj, calls


# -- the marker ------------------------------------------------------------


def test_mapping_writes_the_marker_only_after_every_check_passes(
    tmp_path, monkeypatch
):
    snap, calls = _snap(tmp_path, monkeypatch)
    bam = tmp_path / "out.bam"
    snap.map_reads("l.fq", "r.fq", output=str(bam))

    marker = align_done_path(bam)
    assert marker.exists()
    body = marker.read_text()
    # The command is what lets a reused BAM be matched to the settings that
    # produced it, which is the whole reason the marker beats inferring
    # completion from the file itself.
    assert "snap-aligner paired" in body
    assert "fragments: 500" in body


def test_no_marker_when_snap_exits_zero_having_written_nothing(
    tmp_path, monkeypatch
):
    """The marker must not outlive the checks it is meant to stand for."""
    snap, _ = _snap(tmp_path, monkeypatch, on_run=lambda out: None)
    bam = tmp_path / "out.bam"
    with pytest.raises(Exception, match="produced no alignments"):
        snap.map_reads("l.fq", "r.fq", output=str(bam))
    assert not align_done_path(bam).exists()


def test_alignment_is_done_needs_both_the_marker_and_the_eof(tmp_path):
    bam = tmp_path / "a.bam"
    _complete(bam)
    assert bam_is_complete(bam)
    assert not alignment_is_done(bam), "EOF alone is not the marker"
    align_done_path(bam).write_text("finished\n")
    assert alignment_is_done(bam)

    truncated = tmp_path / "b.bam"
    _partial(truncated)
    align_done_path(truncated).write_text("finished\n")
    assert not alignment_is_done(truncated), "marker alone is not the EOF"


# -- reuse -----------------------------------------------------------------


def test_a_finished_bam_is_reused_rather_than_remapped(tmp_path, monkeypatch):
    snap, calls = _snap(tmp_path, monkeypatch)
    bam = tmp_path / "out.bam"
    _complete(bam)
    align_done_path(bam).write_text("finished\n")
    snap._load_read_count = lambda *a, **k: None

    snap.map_reads("l.fq", "r.fq", output=str(bam))
    assert calls == [], "snap was re-run against a BAM that was already done"


def test_a_bam_predating_the_marker_is_still_reused(tmp_path, monkeypatch):
    """Nobody's existing output directory is invalidated by the marker."""
    snap, calls = _snap(tmp_path, monkeypatch)
    bam = tmp_path / "out.bam"
    _complete(bam)
    snap._load_read_count = lambda *a, **k: None

    snap.map_reads("l.fq", "r.fq", output=str(bam))
    assert calls == [], "a complete BAM was remapped for want of a marker"
    assert align_done_path(bam).exists(), "the marker should be backfilled"


# -- the partial BAM -------------------------------------------------------


def test_a_partial_bam_is_moved_aside_not_overwritten(tmp_path, monkeypatch):
    """The expensive one: this file used to be written straight over.

    It cannot be used for metrics, but it is the only record of what the
    aligner did before it died -- which is exactly what is wanted when the
    crash is amplab/snap#171.
    """
    snap, calls = _snap(tmp_path, monkeypatch)
    bam = tmp_path / "out.bam"
    _partial(bam, size=500)
    original = bam.read_bytes()

    snap.map_reads("l.fq", "r.fq", output=str(bam))

    kept = Path(str(bam) + ".partial")
    assert kept.exists(), "the partial BAM was destroyed"
    assert kept.read_bytes() == original, "the partial BAM was not kept intact"
    assert calls, "snap should have mapped again"
    assert bam_is_complete(bam), "the new BAM should be complete"


def test_only_one_partial_is_kept(tmp_path, monkeypatch):
    """Bounded at one extra file, so retries cannot fill the disk."""
    bam = tmp_path / "out.bam"
    for _ in range(3):
        snap, _ = _snap(tmp_path, monkeypatch)
        _partial(bam)
        snap.map_reads("l.fq", "r.fq", output=str(bam))

    partials = list(tmp_path.glob("*.partial*"))
    assert len(partials) == 1, partials


# -- the run's own lifecycle ----------------------------------------------


def _main(monkeypatch, tmp_path, fail=False):
    """Drive cli.main with the analysis itself stubbed out.

    What is under test is when main() deletes, not what analyse_assembly
    computes, so the assembly step is replaced by something that leaves the
    same files behind.
    """
    import pytransrate.cli as cli

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    out = tmp_path / "out"

    handed = {}

    def fake_analyse(assembly_path, args, result_dir, defer_delete=None):
        result_dir.mkdir(parents=True, exist_ok=True)
        bam = result_dir / "reads.bam"
        _complete(bam)
        align_done_path(bam).write_text("finished\n")
        salmon = result_dir / "salmon"
        salmon.mkdir(exist_ok=True)
        (salmon / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        handed["defer_delete"] = defer_delete
        if defer_delete is not None:
            defer_delete.append(str(bam))
        if fail:
            raise CommandError("the step blew up")
        return {"assembly": str(assembly_path)}

    monkeypatch.setattr(cli, "analyse_assembly", fake_analyse)
    monkeypatch.setattr(cli, "check_arguments", lambda args: [fasta])
    monkeypatch.setattr(cli, "write_assemblies_csv", lambda *a, **k: None)
    code = cli.main(["-a", str(fasta), "-o", str(out), "--no-banner"])
    return code, out / "reads.bam", handed


def test_a_failed_run_keeps_its_bam(tmp_path, monkeypatch):
    """The regression this file exists for.

    The BAM used to be deleted at the end of each assembly's analysis, before
    the run had written assemblies.csv -- so a run that died later had already
    thrown away the hours a rerun needed.
    """
    code, bam, handed = _main(monkeypatch, tmp_path, fail=True)
    assert code == 1
    assert bam.exists(), "a failed run deleted the BAM"
    assert (bam.parent / "salmon" / "quant.sf").exists()
    # main() must hand the analysis somewhere to defer deletion to, which is
    # what stops the BAM being removed at the end of one assembly instead of
    # at the end of the run.  Without this the stub above would be testing
    # only itself: the eager delete this regression is about lives inside
    # analyse_assembly, which the stub replaces.
    assert isinstance(handed["defer_delete"], list)


def test_a_successful_run_still_removes_the_bam(tmp_path, monkeypatch):
    """Deferring the delete must not turn into never deleting."""
    code, bam, _ = _main(monkeypatch, tmp_path, fail=False)
    assert code == 0
    assert not bam.exists()
    assert not align_done_path(bam).exists(), "the marker outlived its BAM"


def test_keep_bam_is_unaffected(tmp_path, monkeypatch):
    import pytransrate.cli as cli

    fasta = tmp_path / "a.fa"
    fasta.write_text(">c\nACGT\n")
    out = tmp_path / "out"
    seen = {}

    def fake_analyse(assembly_path, args, result_dir, defer_delete=None):
        result_dir.mkdir(parents=True, exist_ok=True)
        bam = result_dir / "reads.bam"
        _complete(bam)
        seen["keep"] = args.keep_bam
        if not args.keep_bam and defer_delete is not None:
            defer_delete.append(str(bam))
        return {"assembly": str(assembly_path)}

    monkeypatch.setattr(cli, "analyse_assembly", fake_analyse)
    monkeypatch.setattr(cli, "check_arguments", lambda args: [fasta])
    monkeypatch.setattr(cli, "write_assemblies_csv", lambda *a, **k: None)
    cli.main(["-a", str(fasta), "-o", str(out), "--no-banner", "--keep-bam"])
    assert seen["keep"] is True
    assert (out / "reads.bam").exists()


# -- the inventory ---------------------------------------------------------


def test_a_failure_names_what_it_kept(tmp_path, monkeypatch, caplog):
    """The log has to be able to answer "do I have to map again?"."""
    with caplog.at_level(logging.ERROR, logger="pytransrate"):
        _main(monkeypatch, tmp_path, fail=True)
    text = caplog.text
    assert "reads.bam" in text
    assert "quant.sf" in text
    assert "[complete]" in text


def test_the_inventory_marks_a_partial_bam_as_such(tmp_path, caplog):
    from pytransrate.cli import _report_what_survives

    _partial(tmp_path / "half.bam")
    with caplog.at_level(logging.ERROR, logger="pytransrate"):
        _report_what_survives(tmp_path)
    assert "[partial]" in caplog.text


def test_the_inventory_says_so_when_there_is_nothing(tmp_path, caplog):
    from pytransrate.cli import _report_what_survives

    with caplog.at_level(logging.ERROR, logger="pytransrate"):
        _report_what_survives(tmp_path)
    assert "nothing reusable" in caplog.text
