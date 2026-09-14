"""Tests for gzipped input: detection, reading, and the plain copy.

The end of the chain -- gzipped reads through the real snap-aligner -- is in
tests/test_pipeline.py, since only the binary can answer for that.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from pytransrate.assembly import Assembly, parse_fasta
from pytransrate.compression import (
    is_gzip,
    open_binary,
    open_text,
    plain_copy,
    plain_path,
    strip_gzip_suffix,
)

FASTA = ">c1 a description\nACGTACGT\n>c2\nTTTTGGGG\n"


def _gzipped(tmp_path, text, name="a.fa.gz"):
    path = tmp_path / name
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


# ---------------------------------------------------------------------------
# detection
# ---------------------------------------------------------------------------


def test_gzip_is_detected_by_content_not_by_name(tmp_path):
    """A gzip stream is one whatever it is called."""
    path = _gzipped(tmp_path, FASTA, name="misnamed.fa")
    assert is_gzip(path)
    assert list(parse_fasta(path)) == [("c1", "ACGTACGT"), ("c2", "TTTTGGGG")]


def test_a_plain_file_named_gz_is_not_treated_as_gzip(tmp_path):
    """The other half of the same rule, and the one an extension check fails."""
    path = tmp_path / "plain.fa.gz"
    path.write_text(FASTA)
    assert not is_gzip(path)
    assert list(parse_fasta(path)) == [("c1", "ACGTACGT"), ("c2", "TTTTGGGG")]


def test_unreadable_paths_are_not_gzip(tmp_path):
    """The caller opens it properly a moment later and reports it in context."""
    assert not is_gzip(tmp_path / "missing.fa")
    assert not is_gzip(tmp_path)


def test_short_files_are_not_gzip(tmp_path):
    """Fewer than two bytes to compare against the magic number."""
    path = tmp_path / "tiny.fa"
    path.write_bytes(b"\x1f")
    assert not is_gzip(path)


def test_bgzf_reads_as_gzip(tmp_path):
    """bgzip output is a valid gzip stream, so it needs no special case."""
    pysam = pytest.importorskip("pysam")
    plain = tmp_path / "a.fa"
    plain.write_text(FASTA)
    path = tmp_path / "a.fa.gz"
    pysam.tabix_compress(str(plain), str(path), force=True)
    assert is_gzip(path)
    assert list(parse_fasta(path)) == [("c1", "ACGTACGT"), ("c2", "TTTTGGGG")]


# ---------------------------------------------------------------------------
# openers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("compressed", [False, True])
def test_open_text_and_open_binary_read_either_form(tmp_path, compressed):
    if compressed:
        path = _gzipped(tmp_path, FASTA)
    else:
        path = tmp_path / "a.fa"
        path.write_text(FASTA)

    with open_text(path) as handle:
        assert handle.read() == FASTA
    with open_binary(path) as handle:
        assert handle.read() == FASTA.encode()


# ---------------------------------------------------------------------------
# naming
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, expected",
    [
        ("a.fa.gz", "a.fa"),
        ("a.fa.gzip", "a.fa"),
        ("a.fa.GZ", "a.fa"),
        ("a.fa", "a.fa"),
        ("a.gz.fa", "a.gz.fa"),
    ],
)
def test_strip_gzip_suffix(name, expected):
    assert strip_gzip_suffix(f"/tmp/{name}").name == expected


def test_strip_gzip_suffix_keeps_the_directory():
    assert strip_gzip_suffix("/a/b/c.fa.gz") == Path("/a/b/c.fa")


# ---------------------------------------------------------------------------
# the copy the external tools are given
# ---------------------------------------------------------------------------


def test_plain_copy_writes_the_decompressed_file(tmp_path):
    source = _gzipped(tmp_path, FASTA)
    out = tmp_path / "out"
    out.mkdir()

    copy = plain_copy(source, out)

    assert copy == out / "a.fa"
    assert copy.read_text() == FASTA


def test_plain_copy_refuses_to_overwrite(tmp_path):
    source = _gzipped(tmp_path, FASTA)
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.fa").write_text("someone else's file")

    with pytest.raises(FileExistsError, match="already exists"):
        plain_copy(source, out)
    assert (out / "a.fa").read_text() == "someone else's file"


def test_plain_path_leaves_uncompressed_input_alone(tmp_path):
    """The common case must not copy the assembly at all."""
    source = tmp_path / "a.fa"
    source.write_text(FASTA)
    out = tmp_path / "out"
    out.mkdir()

    with plain_path(source, out) as path:
        assert path == str(source)
    assert list(out.iterdir()) == []


def test_plain_path_decompresses_then_cleans_up(tmp_path):
    source = _gzipped(tmp_path, FASTA)
    out = tmp_path / "out"
    out.mkdir()

    with plain_path(source, out) as path:
        assert Path(path).read_text() == FASTA
        inside = Path(path)
    assert not inside.exists()


def test_plain_path_cleans_up_when_the_run_fails(tmp_path):
    """snap and salmon do fail, and the copy must not outlive them."""
    source = _gzipped(tmp_path, FASTA)
    out = tmp_path / "out"
    out.mkdir()

    with pytest.raises(RuntimeError):
        with plain_path(source, out) as path:
            inside = Path(path)
            raise RuntimeError("snap failed")
    assert not inside.exists()


# ---------------------------------------------------------------------------
# through the loaders
# ---------------------------------------------------------------------------


def test_assembly_loads_gzipped_fasta(tmp_path):
    plain = tmp_path / "a.fa"
    plain.write_text(FASTA)
    gz = _gzipped(tmp_path, FASTA)

    assert Assembly(gz).basic_stats() == Assembly(plain).basic_stats()


def test_read_count_fallback_counts_gzipped_fastq(tmp_path):
    """Counting the lines of a gzip stream would count compressed bytes."""
    from pytransrate.mapper import Snap

    reads = "".join(
        f"@r{i}\nACGT\n+\nIIII\n" for i in range(250)
    )
    path = _gzipped(tmp_path, reads, name="r.fq.gz")

    snap = Snap.__new__(Snap)  # bypass the PATH lookup
    snap._read_count_file = None
    snap._load_read_count(str(path))

    assert snap.read_count == 250


# ---------------------------------------------------------------------------
# a whole sequence-only run
# ---------------------------------------------------------------------------


def test_cli_runs_on_a_gzipped_assembly(tmp_path, monkeypatch):
    """No reads, so no external tool: this is the whole path end to end."""
    from pytransrate.cli import main

    contigs = "".join(f">c{i}\n{'ACGT' * 100}\n" for i in range(3))
    gz = _gzipped(tmp_path, contigs, name="asm.fa.gz")
    out = tmp_path / "results"

    monkeypatch.chdir(tmp_path)
    assert main(["-a", str(gz), "-o", str(out), "--no-banner"]) == 0

    rows = list(csv.reader(open(out / "contigs.csv")))
    assert [row[0] for row in rows[1:]] == ["c0", "c1", "c2"]

    # The run names itself after the assembly, and .gz is not part of that.
    assemblies = list(csv.reader(open(out / "assemblies.csv")))
    name = assemblies[0].index("assembly")
    assert assemblies[1][name] == str(gz)
    assert not list(out.glob("asm.fa"))


def test_cli_result_directory_drops_the_gzip_suffix(tmp_path, monkeypatch):
    """Two assemblies each get a subdirectory named after themselves."""
    from pytransrate.cli import main

    contigs = ">c1\n" + "ACGT" * 100 + "\n"
    one = _gzipped(tmp_path, contigs, name="one.fa.gz")
    two = tmp_path / "two.fa"
    two.write_text(contigs)
    out = tmp_path / "results"

    monkeypatch.chdir(tmp_path)
    assert main(["-a", f"{one},{two}", "-o", str(out), "--no-banner"]) == 0

    assert (out / "one" / "contigs.csv").exists()
    assert (out / "two" / "contigs.csv").exists()
