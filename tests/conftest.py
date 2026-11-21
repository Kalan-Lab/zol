"""Shared pytest configuration and fixtures for the zol test suite."""
import os
import random
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"

# Allow running the tests from a source checkout without installing zol.
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
os.environ["PYTHONPATH"] = os.pathsep.join(
    [str(SRC_DIR)] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _mutate(seq, rng, rate=0.1):
    return "".join(rng.choice(AMINO_ACIDS) if rng.random() < rate else c for c in seq)


def _set_text_attr(obj, attr, value):
    """Set a name/accession on a pyhmmer object (bytes in some versions, str in others)."""
    try:
        setattr(obj, attr, value.encode())
    except TypeError:
        setattr(obj, attr, value)


@pytest.fixture(scope="session")
def hmm_test_data(tmp_path_factory):
    """
    Build a small synthetic HMM database (3 families, FAM0-FAM2, with accessions and trusted / gathering cutoffs)
    and a protein FASTA file with:
      - p1: FAM0 domain (positions 1-120), 4 residue linker, FAM1 domain (positions 125-244)
      - p2: FAM2 domain (positions 1-120)
      - p3: random sequence (no hit expected)
    """
    import pyhmmer

    rng = random.Random(1)
    tmp = tmp_path_factory.mktemp("hmm")
    families = {f"FAM{i}": "".join(rng.choice(AMINO_ACIDS) for _ in range(120)) for i in range(3)}

    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm_file = tmp / "db.hmm"
    with open(hmm_file, "wb") as handle:
        for name, seq in families.items():
            seqs = [
                pyhmmer.easel.TextSequence(name=f"{name}_{j}".encode(), sequence=_mutate(seq, rng))
                for j in range(8)
            ]
            msa = pyhmmer.easel.TextMSA(name=name.encode(), sequences=seqs).digitize(alphabet)
            hmm, _, _ = builder.build_msa(msa, background)
            _set_text_attr(hmm, "accession", name + ".1")
            hmm.cutoffs.trusted = (5.0, 5.0)
            hmm.cutoffs.gathering = (5.0, 5.0)
            hmm.write(handle)

    protein_faa = tmp / "proteins.faa"
    with open(protein_faa, "w") as handle:
        handle.write(">p1\n" + _mutate(families["FAM0"], rng) + "GGGG" + _mutate(families["FAM1"], rng) + "\n")
        handle.write(">p2\n" + _mutate(families["FAM2"], rng) + "\n")
        handle.write(">p3\n" + "".join(rng.choice(AMINO_ACIDS) for _ in range(100)) + "\n")

    empty_faa = tmp / "empty.faa"
    empty_faa.write_text("")

    return {"dir": tmp, "hmm_file": str(hmm_file), "protein_faa": str(protein_faa), "z": len(families)}
