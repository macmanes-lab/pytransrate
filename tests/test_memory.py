"""Tests for the memory budget that caps the assignment workers.

The thing being protected is a nine-hour run: mapping and quantifying finish,
the assignment step sizes its accumulators by the assembly times --threads,
and the OOM killer takes the lot.  See MEMORY_BUDGET in read_metrics.
"""

from __future__ import annotations

import logging

import pytest

from pytransrate import memory
from pytransrate.memory import available_bytes, format_bytes, parse_size
from pytransrate.read_metrics import (
    _MEMORY_HEADROOM,
    _WORKER_OVERHEAD_BYTES,
    _memory_capped_workers,
    _shared_bytes_per_worker,
)

GB = 10 ** 9


# -- parsing --------------------------------------------------------------


@pytest.mark.parametrize(
    "text,expected",
    [
        ("670", 670 * GB),          # a bare number is GB, as schedulers mean
        ("670G", 670 * GB),
        ("670g", 670 * GB),
        ("670GB", 670 * GB),
        ("512M", 512 * 10 ** 6),
        ("1.5T", int(1.5 * 10 ** 12)),
        ("  64G  ", 64 * GB),
        # Binary only when it is asked for, so the log reads back the figure
        # that was typed.
        ("670Gi", 670 * (1 << 30)),
        ("670GiB", 670 * (1 << 30)),
        ("512Mi", 512 * (1 << 20)),
    ],
)
def test_parse_size(text, expected):
    assert parse_size(text) == expected


def test_a_budget_reads_back_as_the_figure_that_was_typed():
    """--mem 670 must not turn into "719.4 GB" in the warning it produces."""
    assert format_bytes(parse_size("670")) == "670.0 GB"


@pytest.mark.parametrize("text", ["", "lots", "-5G", "0", "5X", "G"])
def test_parse_size_rejects_nonsense(text):
    with pytest.raises(ValueError):
        parse_size(text)


def test_format_bytes():
    assert format_bytes(927_700_000_000) == "927.7 GB"


# -- detection ------------------------------------------------------------


def test_available_bytes_takes_the_smallest_limit(monkeypatch):
    """Every limit is real: the node's does not excuse the cgroup's."""
    monkeypatch.setattr(memory, "_cgroup_available", lambda *a: 100 * GB)
    monkeypatch.setattr(memory, "_slurm_allocation", lambda: 40 * GB)
    monkeypatch.setattr(memory, "_meminfo_available", lambda: 900 * GB)
    assert available_bytes() == 40 * GB


def test_available_bytes_is_none_when_nothing_answers(monkeypatch):
    monkeypatch.setattr(memory, "_cgroup_available", lambda *a: None)
    monkeypatch.setattr(memory, "_slurm_allocation", lambda: None)
    monkeypatch.setattr(memory, "_meminfo_available", lambda: None)
    monkeypatch.setattr(memory, "_sysconf_available", lambda: None)
    assert available_bytes() is None


def test_cgroup_limit_read(tmp_path):
    limit = tmp_path / "memory.max"
    usage = tmp_path / "memory.current"
    limit.write_text("64424509440\n")   # 60 GiB, as cgroups count
    usage.write_text("10737418240\n")   # 10 GiB
    assert memory._cgroup_available(limit, usage) == 50 * (1 << 30)


def test_uncapped_cgroup_is_not_a_limit(tmp_path):
    limit = tmp_path / "memory.max"
    usage = tmp_path / "memory.current"
    usage.write_text("0")
    limit.write_text("max")             # cgroup v2
    assert memory._cgroup_available(limit, usage) is None
    limit.write_text(str(2 ** 63 - 4096))  # cgroup v1
    assert memory._cgroup_available(limit, usage) is None


def test_missing_cgroup_files_are_not_a_limit(tmp_path):
    assert memory._cgroup_available(tmp_path / "nope", tmp_path / "nor") is None


def test_slurm_allocation(monkeypatch):
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "256000")
    assert memory._slurm_allocation() == 256000 * (1 << 20)

    monkeypatch.delenv("SLURM_MEM_PER_NODE")
    monkeypatch.setenv("SLURM_MEM_PER_CPU", "4000")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "40")
    assert memory._slurm_allocation() == 4000 * 40 * (1 << 20)


def test_no_slurm_no_allocation(monkeypatch):
    for name in ("SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU",
                 "SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK"):
        monkeypatch.delenv(name, raising=False)
    assert memory._slurm_allocation() is None


# -- capping --------------------------------------------------------------


def test_bytes_per_worker_counts_every_base():
    """4 bytes a base, plus one per contig for the difference array's tail."""
    lengths = [1000, 2000, 500]
    size = _shared_bytes_per_worker(lengths, len(lengths))
    assert size > (sum(lengths) + len(lengths)) * 4


def test_a_budget_that_fits_changes_nothing():
    assert _memory_capped_workers(40, 1 * GB, budget=1000 * GB) == 40


def test_the_reported_failure_is_capped(caplog):
    """40 workers x 23 GB = 928 GB, which is what the OOM killer saw."""
    with caplog.at_level(logging.WARNING, logger="pytransrate"):
        workers = _memory_capped_workers(40, 23 * GB, budget=500 * GB)
    assert 1 < workers < 40
    # It must fit, with the parent's own copy alongside the workers'.
    assert (workers + 1) * 23 * GB <= 500 * GB * _MEMORY_HEADROOM + (
        _WORKER_OVERHEAD_BYTES
    )
    assert "--max-memory" in caplog.text


def test_a_budget_too_small_for_one_worker_falls_back_to_serial(caplog):
    with caplog.at_level(logging.WARNING, logger="pytransrate"):
        assert _memory_capped_workers(40, 200 * GB, budget=100 * GB) == 1
    assert "one process" in caplog.text


def test_an_unknown_budget_does_not_cap(monkeypatch, caplog):
    """A guess that wastes the machine is worse than no guess."""
    monkeypatch.setattr(
        "pytransrate.read_metrics.available_budget", lambda: None
    )
    assert _memory_capped_workers(40, 23 * GB) == 40


def test_serial_is_left_alone():
    assert _memory_capped_workers(1, 500 * GB, budget=1) == 1


def test_the_budget_is_consulted_when_none_is_given(monkeypatch):
    from pytransrate.memory import Budget

    monkeypatch.setattr(
        "pytransrate.read_metrics.available_budget",
        lambda: Budget(500 * GB, "cgroup limit"),
    )
    assert _memory_capped_workers(40, 23 * GB) < 40


def test_the_warning_says_where_the_figure_came_from(caplog):
    """773.1 GB is worth trusting if it is the allocation and worth
    overriding if it is the node's free memory; the number alone says
    neither."""
    from pytransrate.memory import Budget

    with caplog.at_level(logging.WARNING, logger="pytransrate"):
        _memory_capped_workers(
            40, 23 * GB, budget=Budget(720 * (1 << 30), "Slurm allocation")
        )
    assert "Slurm allocation is 773.1 GB" in caplog.text


def test_a_plain_integer_budget_is_the_flag(caplog):
    with caplog.at_level(logging.WARNING, logger="pytransrate"):
        _memory_capped_workers(40, 23 * GB, budget=670 * GB)
    assert "--max-memory setting is 670.0 GB" in caplog.text


# -- the notes ------------------------------------------------------------


def test_usage_quotes_the_warning_this_actually_prints():
    """USAGE.md quotes the capping warning verbatim, and a quote that has
    drifted from the code is worse than no quote: it is the line someone
    greps their log for.  This pins the two together."""
    import re
    from pathlib import Path

    from pytransrate.memory import Budget

    usage = Path(__file__).resolve().parents[1] / "USAGE.md"
    if not usage.exists():        # an installed copy without the docs
        pytest.skip("USAGE.md is not part of this checkout")

    quoted = re.search(
        r"\[ WARN\] (40 processes would need.*?)\n```", usage.read_text(), re.S
    )
    assert quoted, "USAGE.md no longer quotes the capping warning"

    records: list[str] = []

    class _Grab(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    log = logging.getLogger("pytransrate")
    handler = _Grab()
    log.addHandler(handler)
    try:
        _memory_capped_workers(
            40,
            int(927.7e9 / 40),
            budget=Budget(720 * (1 << 30), "Slurm allocation"),
        )
    finally:
        log.removeHandler(handler)

    assert " ".join(quoted.group(1).split()) == records[-1]
