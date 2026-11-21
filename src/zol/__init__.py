"""
ZOL: A comprehensive bioinformatics toolkit for gene cluster analysis.

This package provides tools for detecting and analyzing gene clusters, 
co-located sets of genes, e.g. operons, BGCs, etc.
"""

__author__ = "Rauf Salamzade"
__email__ = "salamzader@gmail.com"
__license__ = "BSD-3-Clause"

# Package metadata
PACKAGE_NAME = "zol"
PACKAGE_DESCRIPTION = "Gene cluster analysis toolkit"
PACKAGE_URL = "https://github.com/Kalan-Lab/zol"

import warnings

# Biopython warns when writing GenBank records whose names (e.g. long scaffold IDs) do not fit the standard LOCUS
# line. The longer LOCUS line is still parsed correctly by Biopython, so the warning is muted for all zol programs
# (this module is imported first by every zol program, script and multiprocessing worker).
warnings.filterwarnings("ignore", message="Increasing length of locus line")


def _read_version() -> str:
    """
    Determine the zol version. pyproject.toml is the single source of truth: when running from a source checkout
    (including editable installs) the version is read directly from it, otherwise the version recorded in the
    installed package metadata (which is generated from pyproject.toml at build time) is used.
    """
    import re
    from pathlib import Path

    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    try:
        in_project = False
        for line in pyproject.read_text().splitlines():
            stripped = line.strip()
            if stripped.startswith("["):
                in_project = stripped == "[project]"
                continue
            if in_project:
                match = re.match(r"""version\s*=\s*["']([^"']+)["']""", stripped)
                if match:
                    return match.group(1)
    except OSError:
        pass

    try:
        from importlib.metadata import PackageNotFoundError, version
        return version(PACKAGE_NAME)
    except PackageNotFoundError:
        return "NA"


__version__ = _read_version()
