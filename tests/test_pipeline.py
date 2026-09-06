"""End-to-end run against the real snap-aligner and salmon.

Skipped unless both are on PATH. To run it:

    micromamba create -p ./env -c conda-forge -c bioconda \\
        snap-aligner=2.0.5 salmon=2.7.0
    PATH=$PWD/env/bin:$PATH pytest tests/test_pipeline.py

This is the test that would have caught the salmon 0.8.2 -> 2.x breakage
(``--useErrorModel`` is now a hard error, ``--sampleOut`` a no-op) and the
snap-aligner soft-clipping change. None of it is visible from unit tests.
"""

from __future__ import annotations

import csv
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("snap-aligner") is None or shutil.which("salmon") is None,
    reason="needs snap-aligner and salmon on PATH",
)

READ_LEN = 100
FRAG_LEN = 250
COMP = str.maketrans("ACGT", "TGCA")


def _rc(seq):
    return seq.translate(COMP)[::-1]


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """A small transcriptome with a planted near-duplicate and a dead contig."""
    rng = random.Random(20260906)
    work = tmp_path_factory.mktemp("pipeline")

    transcripts = {
        f"tx{i}": "".join(rng.choice("ACGT") for _ in range(rng.randint(600, 1500)))
        for i in range(12)
    }
    # tx11 is tx0 with a few substitutions: forces genuine multi-mapping.
    copy = list(transcripts["tx0"])
    for _ in range(len(copy) // 200):
        copy[rng.randrange(len(copy))] = rng.choice("ACGT")
    transcripts["tx11"] = "".join(copy)

    fasta = work / "assembly.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in transcripts.items()))

    r1_path, r2_path = work / "reads.1.fq", work / "reads.2.fq"
    count = 0
    with open(r1_path, "w") as r1, open(r2_path, "w") as r2:
        for i, (name, seq) in enumerate(transcripts.items()):
            depth = 1 if name == "tx10" else (30 if i % 2 == 0 else 12)
            for _ in range(depth):
                if len(seq) < FRAG_LEN:
                    continue
                start = rng.randrange(0, len(seq) - FRAG_LEN)
                frag = seq[start:start + FRAG_LEN]
                qual = "I" * READ_LEN
                r1.write(f"@f{count}/1\n{frag[:READ_LEN]}\n+\n{qual}\n")
                r2.write(f"@f{count}/2\n{_rc(frag[-READ_LEN:])}\n+\n{qual}\n")
                count += 1

    return {
        "dir": work,
        "fasta": fasta,
        "left": r1_path,
        "right": r2_path,
        "fragments": count,
        "transcripts": transcripts,
    }


def _env():
    """pytest's `pythonpath` setting does not reach subprocesses."""
    env = dict(os.environ)
    src = str(Path(__file__).resolve().parent.parent / "src")
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{src}{os.pathsep}{existing}" if existing else src
    return env


@pytest.fixture(scope="module")
def run_output(dataset):
    """Run the CLI exactly as ORP invokes it."""
    out = dataset["dir"] / "out"
    result = subprocess.run(
        [
            sys.executable, "-m", "pytransrate.cli",
            "-o", str(out),
            "-t", "2",
            "-a", str(dataset["fasta"]),
            "--left", str(dataset["left"]),
            "--right", str(dataset["right"]),
        ],
        cwd=dataset["dir"],
        capture_output=True,
        text=True,
        env=_env(),
    )
    if result.returncode != 0:
        pytest.fail(f"transrate failed:\n{result.stdout}\n{result.stderr}")
    return out


def test_produces_the_two_csvs(run_output):
    assert (run_output / "assemblies.csv").exists()
    assert (run_output / "contigs.csv").exists()


def test_assemblies_csv_is_readable_the_way_orp_reads_it(run_output):
    """oyster.py: rows[1][36] and rows[1][37]."""
    rows = list(csv.reader(open(run_output / "assemblies.csv")))
    assert rows[0][36] == "score"
    assert rows[0][37] == "optimal_score"
    assert len(rows[0]) == 40

    score = float(rows[1][36])
    optimal = float(rows[1][37])
    assert 0.0 < score <= 1.0
    assert 0.0 < optimal <= 1.0


def test_contigs_csv_score_column_parses(run_output, dataset):
    """pick_best_contigs.py: float(row[8]), keyed on row[0]."""
    rows = list(csv.reader(open(run_output / "contigs.csv")))
    assert rows[0][8] == "score"

    scores = {row[0]: float(row[8]) for row in rows[1:]}
    assert set(scores) == set(dataset["transcripts"])
    assert all(0.0 < s <= 1.0 for s in scores.values())


def test_dead_contig_scores_at_the_floor(run_output):
    """tx10 was simulated at depth 1 and should bottom out."""
    rows = list(csv.reader(open(run_output / "contigs.csv")))
    scores = {row[0]: float(row[8]) for row in rows[1:]}
    assert scores["tx10"] < 0.1
    assert scores["tx10"] < scores["tx1"]


def test_near_duplicate_loses_to_its_original(run_output):
    """tx11 is a degraded copy of tx0.

    Fragments matching both should be assigned to tx0 on edit distance, so
    tx11 ends up with the poorer score. This is the behaviour salmon 0.8.2's
    postSample.bam used to provide.
    """
    rows = list(csv.reader(open(run_output / "contigs.csv")))
    scores = {row[0]: float(row[8]) for row in rows[1:]}
    assert scores["tx0"] > scores["tx11"]


def test_fragment_count_matches_the_library(run_output, dataset):
    rows = list(csv.reader(open(run_output / "assemblies.csv")))
    fragments = int(rows[1][rows[0].index("fragments")])
    assert fragments == dataset["fragments"]


def test_no_post_sample_bam_is_produced(run_output):
    """salmon 2.x cannot write one; assignment happens in-process."""
    assert not list(run_output.rglob("postSample.bam"))


def test_score_optimisation_curve_is_written(run_output):
    curves = list(run_output.glob("*_score_optimisation.csv"))
    assert curves
    rows = list(csv.reader(open(curves[0])))
    assert rows[0] == ["cutoff", "assembly_score"]


def test_rerun_refuses_to_overwrite(dataset, run_output):
    """The Ruby guarded against clobbering assemblies.csv; so do we."""
    result = subprocess.run(
        [
            sys.executable, "-m", "pytransrate.cli",
            "-o", str(run_output),
            "-a", str(dataset["fasta"]),
        ],
        cwd=dataset["dir"],
        capture_output=True,
        text=True,
        env=_env(),
    )
    assert result.returncode == 1
    assert "would be overwritten" in (result.stdout + result.stderr)
