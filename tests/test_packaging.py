"""Tests for version handling, packaging metadata and Python version compatibility of all scripts."""
import ast
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _pyproject_version():
    text = (REPO_ROOT / "pyproject.toml").read_text()
    project_section = text.split("[project]", 1)[1].split("\n[", 1)[0]
    return re.search(r'^version\s*=\s*"([^"]+)"', project_section, re.M).group(1)


def test_package_version_matches_pyproject():
    import zol
    assert zol.__version__ == _pyproject_version()


def test_zol_suite_version_matches_pyproject():
    result = subprocess.run(
        [sys.executable, "-m", "zol.zol_suite", "--version"], capture_output=True, text=True, cwd=REPO_ROOT
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"ZOL Suite {_pyproject_version()}"


def test_util_get_version_matches_pyproject():
    from zol import util
    assert util.get_version() == _pyproject_version()


def _python_files():
    files = list((REPO_ROOT / "src" / "zol").rglob("*.py")) + list((REPO_ROOT / "scripts").glob("*.py"))
    files += [p for p in (REPO_ROOT / "bin").iterdir() if p.is_file()]
    files += [REPO_ROOT / "docker" / "ZOL", REPO_ROOT / "setup.py"]
    return sorted(files)


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_scripts_parse_with_running_python(path):
    """Guards against syntax only valid in newer Python versions (e.g. nested same-quote f-strings, 3.12+)."""
    ast.parse(path.read_text(), filename=str(path))


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_no_os_system_calls(path):
    tree = ast.parse(path.read_text())
    calls = [
        node.lineno for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        and node.func.attr == "system" and isinstance(node.func.value, ast.Name) and node.func.value.id == "os"
    ]
    assert calls == [], f"os.system used at lines {calls}; use util.run_cmd_via_subprocess / run_cmd_simple instead"


def test_long_locus_name_warning_is_muted_and_record_round_trips(tmp_path):
    """
    Biopython warns ("Increasing length of locus line ...") when writing GenBanks for scaffolds with long names.
    Importing zol mutes it, including in spawned multiprocessing workers; the long LOCUS line still parses.
    Run with -W error so any remaining warning fails the test.
    """
    import subprocess
    import sys
    import textwrap

    src_dir = Path(__file__).resolve().parents[1] / "src"
    script = tmp_path / "write_long_locus.py"
    script.write_text(textwrap.dedent(f"""
        import multiprocessing, sys
        sys.path.insert(0, {str(src_dir)!r})
        import zol  # noqa: F401
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        NAME = "NZ_JAXYZQ010000123_metagenome_scaffold_000042"  # > 16 chars: lengthens LOCUS

        def write_and_read(path):
            rec = SeqRecord(Seq("ATG" * 100), id=NAME, name=NAME, annotations={{"molecule_type": "DNA"}})
            SeqIO.write(rec, path, "genbank")
            back = SeqIO.read(path, "genbank")
            assert back.name == NAME and len(back.seq) == 300, (back.name, len(back.seq))
            return "ok"

        if __name__ == "__main__":
            print(write_and_read({str(tmp_path / "main.gbk")!r}))
            ctx = multiprocessing.get_context("spawn")
            with ctx.Pool(1) as pool:
                print(pool.apply(write_and_read, ({str(tmp_path / "worker.gbk")!r},)))
    """))
    proc = subprocess.run([sys.executable, "-W", "error", str(script)], capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == ["ok", "ok"]
    assert "Increasing length of locus line" not in proc.stderr
