"""Running external commands.

Port of ``lib/transrate/cmd.rb``.  The Ruby built shell strings and handed
them to ``Open3.capture3``; here commands are argument lists run without a
shell, so paths containing spaces or shell metacharacters are safe.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass, field

__all__ = ["CommandError", "CommandResult", "run", "which"]

logger = logging.getLogger("transrate")


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

    def check(self, message: str = "") -> "CommandResult":
        """Raise :class:`CommandError` unless the command succeeded."""
        if self.ok:
            return self
        prefix = message or f"{self.args[0]} failed"
        raise CommandError(
            f"{prefix}\ncommand: {' '.join(self.args)}\n{self.stderr.strip()}"
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


def run(args, cwd=None, env=None, stdout_path=None) -> CommandResult:
    """Run a command and capture its output.

    Args:
        args: program and arguments; every element is coerced to ``str`` so
            callers can pass Paths and ints directly.
        cwd: working directory.
        env: environment overrides, merged over the current environment.
        stdout_path: when given, stdout is written here instead of captured.

    Returns:
        A :class:`CommandResult`; call ``.check()`` to turn failure into an
        exception.
    """
    args = [str(a) for a in args]
    logger.debug("running: %s", " ".join(args))

    full_env = None
    if env:
        import os

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

    proc = subprocess.run(args, cwd=cwd, env=full_env, capture_output=True)
    return CommandResult(
        args=args,
        returncode=proc.returncode,
        stdout=proc.stdout.decode("utf-8", "replace"),
        stderr=proc.stderr.decode("utf-8", "replace"),
    )
