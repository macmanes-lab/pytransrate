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
import gzip
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


def _piped(name):
    """``tx3`` -> ``GADU00000003``, an accession-shaped stand-in.

    Zero-padded so the accessions sort in the same order as the transcripts
    they stand for, which keeps the name-ordered tie-break in assign.py
    (score, then prior, then name) deciding multi-mappers the same way in
    both runs -- tx0 against its near-duplicate tx11, in practice.
    """
    return f"GADU{int(name[2:]):08d}"


def _unpiped(identifier):
    """``ENA|GADU00000003|GADU00000003.1`` -> ``tx3``."""
    return f"tx{int(identifier.split('|')[1][4:])}"


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


@pytest.fixture(scope="module")
def gzipped_run_output(dataset):
    """The same run with every input gzipped.

    This is the only place snap-aligner itself is asked whether it reads
    gzipped FASTQ -- pytransrate passes the read files straight through
    rather than decompressing a library that may be tens of gigabytes, so
    nothing short of the binary can confirm it. The assembly takes the other
    path: snap and salmon are handed a decompressed copy.
    """
    work = dataset["dir"] / "gz"
    work.mkdir(exist_ok=True)

    paths = {}
    for key in ("fasta", "left", "right"):
        source = dataset[key]
        target = work / (source.name + ".gz")
        with open(source, "rb") as handle, gzip.open(target, "wb") as out:
            shutil.copyfileobj(handle, out)
        paths[key] = target

    out = work / "out"
    result = subprocess.run(
        [
            sys.executable, "-m", "pytransrate.cli",
            "-o", str(out),
            "-t", "2",
            "-a", str(paths["fasta"]),
            "--left", str(paths["left"]),
            "--right", str(paths["right"]),
        ],
        cwd=work,
        capture_output=True,
        text=True,
        env=_env(),
    )
    if result.returncode != 0:
        pytest.fail(
            f"transrate failed on gzipped input:\n{result.stdout}\n{result.stderr}"
        )
    return out


@pytest.fixture(scope="module")
def piped_run_output(dataset):
    """The same run with ENA/TSA-style deflines.

    Every transcript is renamed to the shape an ENA or TSA download has --
    ``>ENA|GADU00000003|GADU00000003.1 Gadus morhua mRNA``, two pipes and a
    description. Cutting the identifier at the first ``|``, as the Ruby did,
    named every contig ``ENA``; cutting only at the space is what snap
    writes into the BAM header and what salmon reports, and this is the only
    test that asks the binaries rather than taking that on trust.
    """
    work = dataset["dir"] / "piped"
    work.mkdir(exist_ok=True)

    fasta = work / "assembly.fa"
    fasta.write_text(
        "".join(
            f">ENA|{_piped(n)}|{_piped(n)}.1 Gadus morhua mRNA\n{s}\n"
            for n, s in dataset["transcripts"].items()
        )
    )

    out = work / "out"
    result = subprocess.run(
        [
            sys.executable, "-m", "pytransrate.cli",
            "-o", str(out),
            "-t", "2",
            "-a", str(fasta),
            "--left", str(dataset["left"]),
            "--right", str(dataset["right"]),
        ],
        cwd=work,
        capture_output=True,
        text=True,
        env=_env(),
    )
    if result.returncode != 0:
        pytest.fail(
            f"transrate failed on piped identifiers:\n"
            f"{result.stdout}\n{result.stderr}"
        )
    return out


def test_piped_identifiers_keep_their_pipes(piped_run_output, dataset):
    rows = list(csv.reader(open(piped_run_output / "contigs.csv")))
    names = {row[0] for row in rows[1:]}
    assert names == {
        f"ENA|{_piped(n)}|{_piped(n)}.1" for n in dataset["transcripts"]
    }


def test_piped_identifiers_still_join_to_the_bam(piped_run_output, run_output):
    """The failure this guards is silent: names that match nothing in the
    BAM header produce a full contigs.csv of zeroes rather than an error."""
    def by_name(directory):
        rows = list(csv.reader(open(directory / "contigs.csv")))
        header = rows[0]
        return {
            row[0]: dict(zip(header, row)) for row in rows[1:]
        }

    piped = by_name(piped_run_output)
    plain = by_name(run_output)

    # Every contig was given reads, so a zero here is the join failing.
    assert all(float(row["coverage"]) > 0 for row in piped.values())

    # And the scores are unchanged: only the names differ, and _piped keeps
    # them in an order the name-ordered tie-break decides the same way.
    # Compared the way the gzip test compares them, on score alone -- snap
    # and salmon are both multithreaded and promise no more than that.
    assert {
        _unpiped(name): row["score"] for name, row in piped.items()
    } == {name: row["score"] for name, row in plain.items()}


def test_gzipped_input_scores_identically(run_output, gzipped_run_output):
    """Compression is not allowed to change a single number."""
    def scores(directory):
        rows = list(csv.reader(open(directory / "contigs.csv")))
        return {row[0]: row[8] for row in rows[1:]}

    assert scores(gzipped_run_output) == scores(run_output)


def test_gzipped_assembly_leaves_no_decompressed_copy(gzipped_run_output):
    """The plain copy exists for snap and salmon and no longer."""
    assert not list(gzipped_run_output.rglob("assembly.fa"))


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
