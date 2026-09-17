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

GB = 1 << 30


# -- parsing --------------------------------------------------------------


@pytest.mark.parametrize(
    "text,expected",
    [
        ("200", 200 * GB),          # a bare number is GB, as schedulers mean
        ("200G", 200 * GB),
        ("200g", 200 * GB),
        ("200GB", 200 * GB),
        ("200GiB", 200 * GB),
        ("512M", 512 * (1 << 20)),
        ("1.5T", int(1.5 * (1 << 40))),
        ("  64G  ", 64 * GB),
    ],
)
def test_parse_size(text, expected):
    assert parse_size(text) == expected


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
    limit.write_text("64424509440\n")   # 60 GB
    usage.write_text("10737418240\n")   # 10 GB
    assert memory._cgroup_available(limit, usage) == 50 * GB


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
    monkeypatch.setattr("pytransrate.read_metrics.available_bytes", lambda: None)
    assert _memory_capped_workers(40, 23 * GB) == 40


def test_serial_is_left_alone():
    assert _memory_capped_workers(1, 500 * GB, budget=1) == 1


def test_the_budget_is_consulted_when_none_is_given(monkeypatch):
    monkeypatch.setattr(
        "pytransrate.read_metrics.available_bytes", lambda: 500 * GB
    )
    assert _memory_capped_workers(40, 23 * GB) < 40
