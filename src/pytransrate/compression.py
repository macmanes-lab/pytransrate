"""Reading gzipped inputs.

Assemblies and read files arrive compressed more often than not -- an SRA
download is ``.fastq.gz``, and an assembly kept for any length of time
usually is too -- so both are accepted here in either form.

Compression is detected by **magic number, not by file extension**.  A
gzipped file is not obliged to be named ``.gz``, and a file named ``.gz`` is
not obliged to be gzipped; the first two bytes settle it either way, which is
what samtools and pigz do.  This also means BGZF input works, since a bgzip
file is a valid gzip stream.

Two kinds of consumer need different things:

* Python reads the assembly FASTA and samples the FASTQ itself, so those need
  only an opener that decompresses as it goes -- :func:`open_text`.
* snap-aligner and salmon are handed a *path*, not a file object, so a
  gzipped assembly has to be written out plainly for them first --
  :func:`plain_copy`.  Read files are passed to snap-aligner as they are; it
  reads gzipped FASTQ natively, and decompressing a library that is routinely
  tens of gigabytes would cost more than the whole rest of the run.
"""

from __future__ import annotations

import contextlib
import gzip
import shutil
from pathlib import Path

__all__ = ["GZIP_MAGIC", "is_gzip", "open_text", "open_binary",
           "strip_gzip_suffix", "plain_copy", "plain_path"]

#: The first two bytes of any gzip stream (RFC 1952).
GZIP_MAGIC = b"\x1f\x8b"

#: Suffixes dropped when naming the decompressed form of a file.
_GZIP_SUFFIXES = (".gz", ".gzip")


def is_gzip(path) -> bool:
    """Whether ``path`` is a gzip stream, by its first two bytes.

    Returns ``False`` for anything that cannot be read -- a missing or
    unreadable file is not this function's error to raise, and the caller is
    about to open it properly and report it in context.
    """
    try:
        with open(path, "rb") as handle:
            return handle.read(2) == GZIP_MAGIC
    except OSError:
        return False


def open_text(path):
    """Open ``path`` for reading text, decompressing gzip transparently."""
    if is_gzip(path):
        return gzip.open(path, "rt")
    return open(path)


def open_binary(path):
    """Open ``path`` for reading bytes, decompressing gzip transparently."""
    if is_gzip(path):
        return gzip.open(path, "rb")
    return open(path, "rb")


def strip_gzip_suffix(path):
    """``path`` without a trailing ``.gz``/``.gzip``, as a :class:`Path`.

    Names the file the user would have had uncompressed, which is what the
    output directory and the index should be named after -- ``asm.fa.gz``
    otherwise yields an ``asm.fa`` result directory.
    """
    path = Path(path)
    for suffix in _GZIP_SUFFIXES:
        if path.name.lower().endswith(suffix):
            return path.with_name(path.name[: -len(suffix)])
    return path


def plain_copy(path, directory) -> Path:
    """Decompress ``path`` into ``directory``, returning the new path.

    For the external tools, which take a filename.  The copy is named after
    the source with the gzip suffix removed, and is the caller's to delete.

    Raises:
        FileExistsError: if that name is already taken, rather than
            overwriting something the run did not create.
    """
    target = Path(directory) / strip_gzip_suffix(path).name
    if target.exists():
        raise FileExistsError(
            f"cannot decompress {path} to {target}: the file already exists"
        )
    with gzip.open(path, "rb") as source, open(target, "wb") as handle:
        shutil.copyfileobj(source, handle)
    return target


@contextlib.contextmanager
def plain_path(path, directory):
    """Yield a path to ``path``'s contents that any tool can open.

    Uncompressed input is yielded unchanged, so the common case costs
    nothing.  Gzipped input is decompressed into ``directory`` for the
    duration of the block and removed afterwards -- it is a second copy of
    the input on disk, and nothing outside the block wants it.
    """
    if not is_gzip(path):
        yield str(path)
        return

    target = plain_copy(path, directory)
    try:
        yield str(target)
    finally:
        target.unlink(missing_ok=True)
