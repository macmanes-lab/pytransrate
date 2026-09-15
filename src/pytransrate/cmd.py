"""Running external commands.

Port of ``lib/transrate/cmd.rb``.  The Ruby built shell strings and handed
them to ``Open3.capture3``; here commands are argument lists run without a
shell, so paths containing spaces or shell metacharacters are safe.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["CommandError", "CommandResult", "run", "which"]

logger = logging.getLogger("pytransrate")


class CommandError(Exception):
    """Raised when a required external tool is missing or fails."""


@dataclass
class CommandResult:
    """Outcome of one external command."""

    args: list[str]
    returncode: int
    stdout: str = ""
    stderr: str = field(default="", repr=False)

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def output(self) -> str:
        """Everything the command said, stdout first.

        Callers that only want to search the output for a message should use
        this rather than picking a stream: when ``run`` is given a
        ``log_path`` both streams land in one file and come back together in
        ``stdout``, so ``stderr`` is empty even for a command that failed
        loudly.
        """
        return self.stdout + self.stderr

    def check(self, message: str = "") -> "CommandResult":
        """Raise :class:`CommandError` unless the command succeeded."""
        if self.ok:
            return self
        prefix = message or f"{self.args[0]} failed"
        raise CommandError(
            f"{prefix}\ncommand: {' '.join(self.args)}\n{self.output.strip()}"
        )

    def __str__(self) -> str:
        return " ".join(self.args)


def which(binary: str) -> str:
    """Locate an executable, or raise.

    Replaces the Ruby's ``Cmd.new("which #{bin}")`` idiom.
    """
    path = shutil.which(binary)
    if path is None:
        raise CommandError(f"could not find {binary} in PATH")
    return path


def run(args, cwd=None, env=None, stdout_path=None, log_path=None) -> CommandResult:
    """Run a command and capture its output.

    Args:
        args: program and arguments; every element is coerced to ``str`` so
            callers can pass Paths and ints directly.
        cwd: working directory.
        env: environment overrides, merged over the current environment.
        stdout_path: when given, stdout is written here instead of captured.
        log_path: when given, both streams are appended to this file *as the
            command runs* and returned together in ``stdout``; see
            LIVE_LOGGING.

    Returns:
        A :class:`CommandResult`; call ``.check()`` to turn failure into an
        exception.
    """
    args = [str(a) for a in args]
    logger.debug("running: %s", " ".join(args))

    full_env = None
    if env:
        full_env = {**os.environ, **{k: str(v) for k, v in env.items()}}

    if stdout_path is not None:
        with open(stdout_path, "wb") as handle:
            proc = subprocess.run(
                args, cwd=cwd, env=full_env, stdout=handle, stderr=subprocess.PIPE
            )
        return CommandResult(
            args=args,
            returncode=proc.returncode,
            stdout="",
            stderr=proc.stderr.decode("utf-8", "replace"),
        )

    if log_path is not None:
        return _run_logged(args, log_path, cwd, full_env)

    proc = subprocess.run(args, cwd=cwd, env=full_env, capture_output=True)
    return CommandResult(
        args=args,
        returncode=proc.returncode,
        stdout=proc.stdout.decode("utf-8", "replace"),
        stderr=proc.stderr.decode("utf-8", "replace"),
    )


# ---------------------------------------------------------------------------
# LIVE_LOGGING
#
# capture_output holds everything in pipes until the child exits, so a run
# that is killed rather than returning -- the OOM killer taking snap-aligner
# on a large assembly, a scheduler wall clock, a ^C -- leaves nothing behind
# to explain itself.  That is precisely the run whose log is wanted.
#
# _run_logged hands the child a file descriptor instead, so its output is on
# disk the moment it is written and survives any death the parent doesn't.
# stdout and stderr share one descriptor, so the two interleave in the order
# snap actually produced them.  The cost is that they can no longer be told
# apart afterwards; see CommandResult.output.
# ---------------------------------------------------------------------------


def _run_logged(args, log_path, cwd, full_env) -> CommandResult:
    path = Path(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "ab") as handle:
        handle.seek(0, os.SEEK_END)
        start = handle.tell()
        # Record the command too: a log of four snap runs at different
        # -locationSize values is unreadable without it.
        handle.write(f"$ {' '.join(args)}\n".encode())
        handle.flush()
        proc = subprocess.run(
            args, cwd=cwd, env=full_env, stdout=handle, stderr=handle
        )

    with open(path, "rb") as handle:
        handle.seek(start)
        written = handle.read()

    return CommandResult(
        args=args,
        returncode=proc.returncode,
        stdout=written.decode("utf-8", "replace"),
        stderr="",
    )
