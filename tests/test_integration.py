"""End-to-end smoke test over the offline half of the pipeline.

Everything from a FASTA plus a BAM through to the two CSVs ORP reads.  The
aligner and quantifier are not involved -- this exercises the modules that
can be validated without external binaries, and proves they compose.
"""

from __future__ import annotations

import csv

import numpy as np
import pysam
import pytest

from pytransrate.assembly import Assembly
from pytransrate.bam_metrics import compute_bam_metrics
from pytransrate.output import write_assemblies_csv, write_contigs_csv
from pytransrate.score import ScoreOptimiser

CONTIG_LEN = 600
N_CONTIGS = 6
READ_LEN = 100


@pytest.fixture
def workspace(tmp_path):
    """A small assembly with reads tiling every contig but the last two."""
    rng = np.random.default_rng(20260906)

    seqs = {}
    for i in range(N_CONTIGS):
        seqs[f"c{i}"] = "".join(rng.choice(list("ACGT"), size=CONTIG_LEN))

    fasta = tmp_path / "assembly.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))

    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": n, "LN": CONTIG_LEN} for n in seqs],
    }

    bam_path = tmp_path / "aln.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        pair = 0
        for ref_id in range(N_CONTIGS):
            # Last contig gets no reads at all; second-to-last gets one pair.
            if ref_id == N_CONTIGS - 1:
                continue
            starts = [0] if ref_id == N_CONTIGS - 2 else range(0, 380, 20)
            for start in starts:
                mate_start = start + 120
                for is_read1 in (True, False):
                    read = pysam.AlignedSegment()
                    read.query_name = f"frag{pair}"
                    flag = 1 | (2 if True else 0)
                    if is_read1:
                        flag |= 64 | 32
                        pos, mpos = start, mate_start
                    else:
                        flag |= 128 | 16
                        pos, mpos = mate_start, start
                    read.flag = flag
                    read.reference_id = ref_id
                    read.reference_start = pos
                    read.mapping_quality = 60
                    read.cigarstring = "10S90M" if is_read1 else "100M"
                    read.next_reference_id = ref_id
                    read.next_reference_start = mpos
                    read.query_sequence = "A" * READ_LEN
                    read.query_qualities = pysam.qualitystring_to_array(
                        "I" * READ_LEN
                    )
                    read.set_tag("NM", 2, value_type="i")
                    out.write(read)
                pair += 1

    return {"fasta": fasta, "bam": bam_path, "fragments": pair, "dir": tmp_path}


def test_full_offline_pipeline(workspace):
    assembly = Assembly(workspace["fasta"])
    assert len(assembly) == N_CONTIGS

    contigs = compute_bam_metrics(str(workspace["bam"]), realistic_distance=400)
    by_name = {c.name: c for c in contigs}

    total_good = 0
    for name, contig in assembly:
        metrics = by_name[name]
        contig.set_uncovered_bases(metrics.bases_uncovered)
        contig.p_seq_true = metrics.p_seq_true()
        contig.p_not_segmented = metrics.p_not_segmented
        contig.in_bridges = metrics.bridges
        contig.good = metrics.good
        if metrics.fragments_mapped > 1:
            contig.p_good = metrics.good / metrics.fragments_mapped
        contig.tpm = 10.0
        total_good += metrics.good

    # The well-covered contigs should outscore the barely-covered one, and
    # the unmapped contig should sit at the floor.
    assert assembly["c0"].score > assembly["c4"].score
    assert assembly["c5"].score == pytest.approx(0.01)

    optimiser = ScoreOptimiser(
        assembly=assembly, fragments=workspace["fragments"], good=total_good
    )
    raw = optimiser.raw_score()
    optimal, cutoff = optimiser.optimal_score(
        csv_path=str(workspace["dir"] / "opt.csv")
    )

    assert 0.0 < raw <= 1.0
    # Dropping the dead contigs must beat keeping them.
    assert optimal > raw

    # CUTOFF_BOUNDARY: the dead contigs sit exactly on the returned cutoff
    # (both at the 0.01 floor), and classify uses >=, so they survive even
    # though the optimum was reached by removing them.
    assert cutoff == pytest.approx(0.01)
    assembly.classify_contigs(cutoff)
    assert assembly.good_contigs == N_CONTIGS

    # Anywhere above the floor the split is real.
    assembly.classify_contigs(0.5)
    assert assembly.good_contigs == 4

    contigs_csv = workspace["dir"] / "contigs.csv"
    write_contigs_csv(assembly, str(contigs_csv))

    result = {"assembly": str(workspace["fasta"])}
    result.update(assembly.basic_stats())
    result.update(assembly.contig_metrics())
    result["score"] = raw
    result["optimal_score"] = optimal
    result["cutoff"] = cutoff
    result["weighted"] = optimiser.weighted_score()

    assemblies_csv = workspace["dir"] / "assemblies.csv"
    write_assemblies_csv([result], str(assemblies_csv))

    # Read them back exactly the way ORP does.
    rows = list(csv.reader(open(contigs_csv)))
    assert rows[0][8] == "score"
    scores = {row[0]: float(row[8]) for row in rows[1:]}
    assert len(scores) == N_CONTIGS
    assert scores["c0"] > scores["c5"]

    rows = list(csv.reader(open(assemblies_csv)))
    assert float(rows[1][36]) == pytest.approx(round(raw, 5))
    assert float(rows[1][37]) == pytest.approx(round(optimal, 5))


def test_soft_clips_do_not_leak_into_coverage(workspace):
    """Read 1 carries 10S90M throughout; coverage must still start at POS."""
    contigs = compute_bam_metrics(str(workspace["bam"]), realistic_distance=400)
    coverage = {c.name: np.asarray(c.coverage) for c in contigs}
    # c0's first read pair starts at 0, so base 0 must be covered.
    assert coverage["c0"][0] > 0
    assert coverage["c5"].sum() == 0
