"""Assembly container and assembly-level statistics.

Port of ``lib/transrate/assembly.rb`` and ``lib/transrate/contig_metrics.rb``.

``basic_stats`` reproduces ``basic_bin_stats`` including its several
oddities, which are marked STATS_QUIRK below.  Key order is load-bearing:
ORP reads ``assemblies.csv`` by column index (score at 36, optimal at 37 --
see oyster.py), so the order here is part of the output contract.
"""

from __future__ import annotations

from collections import OrderedDict

from transrate.contig import Contig

__all__ = ["AssemblyError", "Assembly", "BASIC_STATS_KEYS", "CONTIG_METRICS_KEYS"]


class AssemblyError(Exception):
    """Raised for malformed or ambiguous assemblies."""


#: Nx values reported, in the order the Ruby emits them.
_NX = (90, 70, 50, 30, 10)

#: ``basic_stats`` key order.
BASIC_STATS_KEYS = (
    "n_seqs",
    "smallest",
    "largest",
    "n_bases",
    "mean_len",
    "n_under_200",
    "n_over_1k",
    "n_over_10k",
    "n_with_orf",
    "mean_orf_percent",
) + tuple(f"n{n}" for n in _NX)

#: ``ContigMetrics#results`` key order.
CONTIG_METRICS_KEYS = ("gc", "bases_n", "proportion_n")

#: Contigs shorter than this are excluded from the length statistics.
MIN_CONTIG_LENGTH = 200

#: An ORF must exceed this many codons to count toward ``n_with_orf``.
MIN_ORF_CODONS = 149


def parse_fasta(path):
    """Yield ``(name, sequence)``.

    The identifier is everything up to the first whitespace or ``|``,
    matching BioRuby's ``entry_id``, which is what the Ruby used for contig
    names.
    """
    name = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                defline = line[1:].strip()
                token = defline.split()[0] if defline else ""
                name = token.split("|")[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


class Assembly:
    """A transcriptome assembly loaded from FASTA."""

    def __init__(self, path):
        self.file = str(path)
        self.contigs: "OrderedDict[str, Contig]" = OrderedDict()
        self.n_bases = 0
        self._basic_stats = None

        for name, seq in parse_fasta(self.file):
            if len(seq) == 0:
                raise AssemblyError(f"Entry found with no sequence: {name}")
            contig = Contig(name=name, seq=seq)
            if contig.name in self.contigs:
                raise AssemblyError(
                    f"Non-unique fasta identifier found: >{contig.name}\n"
                    "Contig names are taken from before the first | or space. "
                    "If you used Trinity, replace | with _ in the identifiers."
                )
            if "," in contig.name:
                raise AssemblyError("Contig names can't contain commas")
            self.n_bases += contig.length
            self.contigs[contig.name] = contig

    def __len__(self) -> int:
        return len(self.contigs)

    def __getitem__(self, name) -> Contig:
        return self.contigs[name]

    def __contains__(self, name) -> bool:
        return name in self.contigs

    def __iter__(self):
        return iter(self.contigs.items())

    @property
    def size(self) -> int:
        return len(self.contigs)

    def values(self):
        return self.contigs.values()

    # ------------------------------------------------------------------
    # STATS_QUIRK
    #
    # basic_bin_stats has several behaviours worth naming, all preserved:
    #
    #  1. Contigs under 200 bases are skipped for the length accumulation
    #     but still counted in n_seqs, smallest, and n_bases -- so mean_len
    #     divides a >=200bp cumulative length by the *full* contig count and
    #     is not the mean of anything in particular.
    #  2. Nx cutoffs are compared against total n_bases (all contigs) while
    #     the running total only includes >=200bp contigs, so on assemblies
    #     with many short contigs the larger Nx values may never be reached.
    #  3. At most one Nx threshold is satisfied per contig, even when a
    #     single long contig spans several. Unfilled slots are padded with
    #     the longest contig's length.
    #  4. Contigs are walked shortest-first, so the value recorded at the
    #     10% cutoff is N90, at 30% is N70, and so on. The labels are right.
    # ------------------------------------------------------------------

    def basic_stats(self) -> "OrderedDict[str, float]":
        """Length and ORF statistics. See STATS_QUIRK."""
        if self._basic_stats is not None:
            return self._basic_stats

        contigs = sorted(self.contigs.values(), key=lambda c: c.length)
        stats: "OrderedDict[str, float]" = OrderedDict()

        if not contigs:
            for key in BASIC_STATS_KEYS:
                stats[key] = 0
            self._basic_stats = stats
            return stats

        pending = list(_NX)  # popped from the end: 10, 30, 50, 70, 90
        cutoff = pending.pop() / 100.0

        cumulative = 0.0
        results: list[int] = []
        n_under_200 = n_over_1k = n_over_10k = n_with_orf = 0
        orf_length_sum = 0

        for contig in contigs:
            if contig.length < MIN_CONTIG_LENGTH:
                n_under_200 += 1
                continue
            if contig.length > 1_000:
                n_over_1k += 1
            if contig.length > 10_000:
                n_over_10k += 1

            orf_length = contig.orf_length
            orf_length_sum += orf_length
            if orf_length > MIN_ORF_CODONS:
                n_with_orf += 1

            cumulative += contig.length
            if cumulative >= self.n_bases * cutoff:
                results.append(contig.length)
                cutoff = pending.pop() / 100.0 if pending else 1.0

        # Pad unreached thresholds with the longest contig.
        while len(results) < len(_NX):
            results.append(contigs[-1].length)

        mean = cumulative / self.size if self.size else 0.0
        if self.size * mean == 0:
            mean_orf_percent = 0.0
        else:
            mean_orf_percent = 300 * orf_length_sum / (self.size * mean)

        stats["n_seqs"] = len(contigs)
        stats["smallest"] = contigs[0].length
        stats["largest"] = contigs[-1].length
        stats["n_bases"] = self.n_bases
        stats["mean_len"] = mean
        stats["n_under_200"] = n_under_200
        stats["n_over_1k"] = n_over_1k
        stats["n_over_10k"] = n_over_10k
        stats["n_with_orf"] = n_with_orf
        stats["mean_orf_percent"] = mean_orf_percent
        for label, value in zip(_NX, results):
            stats[f"n{label}"] = value

        self._basic_stats = stats
        return stats

    def contig_metrics(self) -> "OrderedDict[str, float]":
        """Assembly-wide base composition. Port of ContigMetrics#results."""
        total = a = c = g = t = bases_n = 0
        for contig in self.contigs.values():
            total += contig.length
            a += contig.bases_a
            c += contig.bases_c
            g += contig.bases_g
            t += contig.bases_t
            bases_n += contig.bases_n

        acgt = a + c + g + t
        out: "OrderedDict[str, float]" = OrderedDict()
        # Note the differing denominators: gc excludes N, proportion_n does not.
        out["gc"] = (g + c) / acgt if acgt else 0.0
        out["bases_n"] = bases_n
        out["proportion_n"] = bases_n / total if total else 0.0
        return out

    def classify_contigs(self, cutoff: float) -> None:
        for contig in self.contigs.values():
            contig.classify(cutoff)

    @property
    def good_contigs(self) -> int:
        return sum(
            1 for c in self.contigs.values() if c.classification == "good"
        )
