"""The startup banner.

Keeps the shape of the Ruby's banner -- same figlet font, same coloured
``░▓▓▓^▓▓▓░`` flanks -- so the lineage is visible, while the wordmark and
version make clear this is not the original.

The banner goes to stderr: it is decoration, not data.  The run report --
every INFO line, including the contig and mapping metric blocks -- goes to
stdout instead, so redirecting stdout captures the numbers and nothing else.
"""

from __future__ import annotations

import os
import shutil
import sys

from pytransrate import __version__

__all__ = ["TAGLINE", "banner", "print_banner"]

#: One-line description of what the tool is for.
TAGLINE = "quality assessment of de-novo transcriptome assemblies"

_ATTRIBUTION = "Python port of transrate, Smith-Unna et al. 2016"

# figlet "standard", the font the Ruby used for its own wordmark.
_WORDMARK = r"""
             _                                 _
 _ __  _   _| |_ _ __ __ _ _ __  ___ _ __ __ _| |_ ___
| '_ \| | | | __| '__/ _` | '_ \/ __| '__/ _` | __/ _ \
| |_) | |_| | |_| | | (_| | | | \__ \ | | (_| | ||  __/
| .__/ \__, |\__|_|  \__,_|_| |_|___/_|  \__,_|\__\___|
|_|    |___/
""".strip("\n")

#: The Ruby's flanking motif, from http://xkcd.com/1179/ by way of snap.rb.
_FLANK = "░▓▓▓^▓▓▓░"

_GREEN, _YELLOW, _RED, _DIM, _RESET = (
    "\033[32m", "\033[33m", "\033[31m", "\033[2m", "\033[0m"
)

#: Below this many columns the full wordmark wraps and looks broken.
_MIN_WIDTH = 72


def _use_colour(stream) -> bool:
    """Colour only for a real terminal that has not opted out."""
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("TERM") == "dumb":
        return False
    return bool(getattr(stream, "isatty", lambda: False)())


def _terminal_width(default: int = 80) -> int:
    try:
        return shutil.get_terminal_size((default, 24)).columns
    except Exception:
        return default


def banner(colour: bool = False, width: int | None = None) -> str:
    """Render the banner.

    Args:
        colour: emit ANSI colour for the flanking motif.
        width: terminal width; detected when omitted. Narrow terminals get a
            single compact line instead of the wordmark.

    Returns:
        The banner text, without a trailing newline. Always carries the
        version.
    """
    if width is None:
        width = _terminal_width()

    version = f"v{__version__}"

    if width < _MIN_WIDTH:
        # Still name the tool and the version -- that is the part that
        # matters in a log or a narrow pane.
        return f"pytransrate {version}\n{TAGLINE}"

    lines = _WORDMARK.split("\n")
    # Indent the unflanked lines so the wordmark stays aligned with the
    # flanked body, as the Ruby's banner did.
    pad = " " * (len(_FLANK) + 1)
    for offset in (0, 1, 5):
        lines[offset] = pad + lines[offset]

    if colour:
        flanks = [_GREEN + _FLANK + _RESET,
                  _YELLOW + _FLANK + _RESET,
                  _RED + _FLANK + _RESET]
        # Flank the three body lines, as the Ruby did.
        for offset, flank in zip(range(2, 5), flanks):
            lines[offset] = f"{flank} {lines[offset]} {flank}"
    else:
        for offset in range(2, 5):
            lines[offset] = f"{_FLANK} {lines[offset]} {_FLANK}"

    footer = f"  {version}  ·  {TAGLINE}"
    attribution = f"  {_ATTRIBUTION}"
    if colour:
        footer = f"  {version}  ·  {TAGLINE}"
        attribution = f"{_DIM}  {_ATTRIBUTION}{_RESET}"

    return "\n".join(["", *lines, "", footer, attribution, ""])


def print_banner(stream=None) -> None:
    """Write the banner to stderr (or ``stream``)."""
    stream = stream or sys.stderr
    print(banner(colour=_use_colour(stream)), file=stream)
