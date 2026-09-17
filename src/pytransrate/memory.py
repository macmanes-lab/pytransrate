"""How much memory this run may actually use.

The read-metrics step sizes its shared accumulators by the assembly and then
multiplies them by ``--threads`` (see MEMORY_BUDGET in
:mod:`~pytransrate.read_metrics`), so on a large assembly the difference
between a run that finishes and a run the OOM killer takes out is a number
nobody can work out in their head.  This module works it out instead.

Nothing here is authoritative -- it is a budget, not a guarantee, and every
source below can be absent.  ``available_bytes`` returns ``None`` when it
cannot tell, and callers are expected to carry on rather than refuse to run.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path

__all__ = ["available_bytes", "format_bytes", "parse_size"]

logger = logging.getLogger("pytransrate")

#: cgroup v2, then v1.  A container or a scheduler-managed job is capped well
#: below what /proc/meminfo reports, because /proc is the host's.
_CGROUP_V2 = (Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current"))
_CGROUP_V1 = (
    Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    Path("/sys/fs/cgroup/memory/memory.usage_in_bytes"),
)

#: Above this, a cgroup limit is "no limit" expressed as a huge number --
#: v1 writes 2**63 rounded down to a page, and anything near it is not a cap.
_NO_LIMIT = 2 ** 62

_SIZE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*([kmgtp]?)i?b?\s*$", re.IGNORECASE)

_UNITS = {"": 1 << 30, "k": 1 << 10, "m": 1 << 20, "g": 1 << 30, "t": 1 << 40,
          "p": 1 << 50}


def parse_size(text) -> int:
    """Bytes from a human-written size: ``200G``, ``512M``, ``1.5T``, ``200``.

    A bare number is gigabytes, because that is the unit every scheduler
    directive and every person asking for memory on a cluster uses.
    """
    match = _SIZE.match(str(text))
    if not match:
        raise ValueError(
            f"cannot read {text!r} as a memory size; try something like 200G"
        )
    value, unit = match.groups()
    size = int(float(value) * _UNITS[unit.lower()])
    if size <= 0:
        raise ValueError(f"memory size must be positive, got {text!r}")
    return size


def format_bytes(size) -> str:
    """``size`` in GB, for a log line."""
    return f"{size / 1e9:.1f} GB"


def _read_int(path):
    """The single integer in ``path``, or None if it is not one."""
    try:
        text = path.read_text().strip()
    except OSError:
        return None
    try:
        return int(text)
    except ValueError:  # cgroup v2 writes "max" for no limit
        return None


def _cgroup_available(limit_path, usage_path):
    """Headroom left in a cgroup, or None if it is not capped by one."""
    limit = _read_int(limit_path)
    if limit is None or limit >= _NO_LIMIT:
        return None
    usage = _read_int(usage_path) or 0
    return max(limit - usage, 0)


def _meminfo_available():
    """``MemAvailable`` -- what can be had without swapping, per the kernel."""
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except (OSError, IndexError, ValueError):
        return None
    return None


def _sysconf_available():
    """Free physical pages.  The fallback where there is no /proc (macOS)."""
    try:
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (AttributeError, ValueError, OSError):
        return None


def _slurm_allocation():
    """What Slurm gave this job, which is what it will enforce.

    A cgroup usually enforces it too, but not every site configures
    ``ConstrainRAMSpace``, and where it does not the job is killed by the
    scheduler against these numbers while /proc/meminfo reports the whole
    node.  Megabytes, as ``sbatch --mem`` takes them.
    """
    per_node = os.environ.get("SLURM_MEM_PER_NODE")
    if per_node:
        try:
            return int(per_node) * (1 << 20)
        except ValueError:
            return None

    per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    cpus = os.environ.get("SLURM_CPUS_ON_NODE") or os.environ.get(
        "SLURM_CPUS_PER_TASK"
    )
    if per_cpu and cpus:
        try:
            return int(per_cpu) * int(cpus) * (1 << 20)
        except ValueError:
            return None
    return None


def available_bytes():
    """Memory this process may reasonably expect to get, or None.

    The smallest of every limit that can be found, because they are all real
    -- a cgroup does not care that the node has more, and the node does not
    care that the cgroup allows more.
    """
    found = {
        "cgroup v2": _cgroup_available(*_CGROUP_V2),
        "cgroup v1": _cgroup_available(*_CGROUP_V1),
        "slurm": _slurm_allocation(),
        "meminfo": _meminfo_available(),
    }
    known = {name: size for name, size in found.items() if size}
    if not known:
        fallback = _sysconf_available()
        if not fallback:
            return None
        known["sysconf"] = fallback

    name = min(known, key=lambda key: known[key])
    logger.debug(
        "memory available: %s (%s); saw %s",
        format_bytes(known[name]),
        name,
        ", ".join(f"{k}={format_bytes(v)}" for k, v in known.items()),
    )
    return known[name]
