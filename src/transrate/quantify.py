"""Transcript quantification with salmon.

Port of ``lib/transrate/salmon.rb``, targeting **salmon 2.7.0** -- the Rust
rewrite -- rather than 0.8.2.

Three changes were forced by the rewrite, all verified against the binary:

1.  ``--useErrorModel`` is a hard error in 2.x ("unexpected argument ...
    tip: a similar argument exists: '--errorModel'").  See ERROR_MODEL.
2.  ``--sampleOut`` / ``--sampleUnaligned`` are accepted but inert: salmon
    logs "accepted but not yet implemented and have no effect" and writes no
    BAM.  ``postSample.bam`` no longer exists, so fragment assignment moved
    in-process to :mod:`transrate.assign`.
3.  ``--libType a`` became ``-l A``.

``quant.sf`` itself is unchanged -- still Name, Length, EffectiveLength,
TPM, NumReads -- so the parser below matches the Ruby's.
"""

from __future__ import annotations

import logging
from pathlib import Path

from transrate.cmd import CommandError, run, which

__all__ = ["SalmonError", "Salmon", "Expression", "load_expression"]

logger = logging.getLogger("transrate")


class SalmonError(CommandError):
    """salmon failed."""


#: quant.sf column count, as the Ruby also asserted.
_QUANT_COLUMNS = 5


# ---------------------------------------------------------------------------
# ERROR_MODEL
#
# salmon 2.6 made deterministic quantification the default and switched
# alignment-mode scoring to the BAM `AS` tag. snap-aligner emits no AS tag --
# only NM -- so on a snap BAM salmon warns:
#
#   "no alignment in this BAM carries an AS tag (910 records read): scoring
#    falls back to equal weights for every placement of a fragment, so
#    multireads are apportioned by effective length alone; pass --errorModel
#    with -t <transcripts.fa> to score the alignments from their own bases"
#
# Apportioning multireads by effective length alone is materially worse on a
# de-novo transcriptome, which is full of near-duplicate contigs -- precisely
# the fragments whose placement decides contig scores. So --errorModel is
# passed by default here. On the synthetic check it resolved 21 equivalence
# classes against 19 without it.
#
# --errorModel requires -t/--targets, which we already pass.
# ---------------------------------------------------------------------------


class Expression(dict):
    """Per-transcript expression: name -> {eff_len, eff_count, tpm}."""


def load_expression(path) -> Expression:
    """Parse ``quant.sf``.

    Raises:
        SalmonError: if the column count is not 5, which historically meant
            a mismatched salmon version.
    """
    expression = Expression()
    with open(path) as handle:
        next(handle, None)  # header
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != _QUANT_COLUMNS:
                raise SalmonError(
                    f"salmon quant.sf should have {_QUANT_COLUMNS} columns "
                    f"but had {len(fields)}. Check your salmon version."
                )
            name, _length, eff_len, tpm, eff_count = fields
            expression[name] = {
                # The Ruby read effective length as an integer; salmon emits
                # it as a float, so truncation is preserved deliberately to
                # keep the coverage arithmetic identical.
                "eff_len": int(float(eff_len)),
                "eff_count": float(eff_count),
                "tpm": float(tpm),
            }
    return expression


class Salmon:
    """Runs ``salmon quant`` in alignment mode."""

    def __init__(self, binary: str | None = None):
        self.binary = binary or which("salmon")

    def version(self) -> str:
        result = run([self.binary, "--version"])
        return result.stdout.strip() or result.stderr.strip()

    def build_command(
        self,
        assembly,
        bam,
        threads: int = 8,
        output_dir: str = "salmon_out",
        lib_type: str = "A",
        error_model: bool = True,
        seq_bias: bool = True,
        gc_bias: bool = True,
    ):
        """Assemble the ``salmon quant`` command.

        Args:
            error_model: pass ``--errorModel``. On by default; see
                ERROR_MODEL for why it matters with snap input.
        """
        args = [
            self.binary, "quant",
            "--alignments", str(bam),
            "--targets", str(assembly),
            "--threads", threads,
            "--output", str(output_dir),
            "--libType", lib_type,
            "--no-version-check",
        ]
        if error_model:
            args.append("--errorModel")
        if seq_bias:
            args.append("--seqBias")
        if gc_bias:
            args.append("--gcBias")
        return args

    def run(
        self,
        assembly,
        bam,
        threads: int = 8,
        output_dir: str = "salmon_out",
        **kwargs,
    ) -> Expression:
        """Quantify, returning parsed expression.

        Reuses an existing ``quant.sf`` if one is already present, matching
        the Ruby's resume behaviour.
        """
        output_dir = Path(output_dir)
        quant_sf = output_dir / "quant.sf"

        if quant_sf.exists():
            logger.info("reusing existing salmon output: %s", quant_sf)
            return load_expression(quant_sf)

        args = self.build_command(
            assembly, bam, threads=threads, output_dir=output_dir, **kwargs
        )
        result = run(args)
        if not result.ok:
            raise SalmonError(f"salmon failed\n{result.stderr}")

        for line in result.stderr.splitlines():
            if "WARN" in line and "AS tag" in line:
                logger.warning("salmon: alignments carry no AS tag")

        if not quant_sf.exists():
            raise SalmonError(f"salmon produced no quant.sf in {output_dir}")

        return load_expression(quant_sf)
