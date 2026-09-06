"""pytransrate: quality assessment of de-novo transcriptome assemblies.

A Python port of transrate (Smith-Unna et al. 2016), targeting current
snap-aligner and salmon. Renamed to keep it distinguishable from the
original Ruby implementation, whose scores it deliberately does not
reproduce -- see the module docstrings for where and why behaviour differs.

Original: https://github.com/blahah/transrate -- see CITATION.md.
"""

__version__ = "2.0.0.dev0"

__all__ = ["__version__"]
