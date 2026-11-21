"""Tests for running external commands and the helpers which replaced os.system calls."""
import gzip
import logging
import os
import subprocess

import pytest

from zol import util


def test_cmd_to_shell_str_quotes_arguments_but_keeps_operators():
    cmd = ["echo", "a b", "it's", "$HOME", ">", "out file.txt", " >> ", "x", ";", "ls", "", "  "]
    assert util.cmd_to_shell_str(cmd) == \
        "echo 'a b' 'it'\"'\"'s' '$HOME' > 'out file.txt' >> x ; ls"


def test_cmd_to_shell_str_passes_strings_through():
    assert util.cmd_to_shell_str("echo hi | cat") == "echo hi | cat"


def test_cmd_to_shell_str_converts_non_strings():
    assert util.cmd_to_shell_str(["head", "-n", 5]) == "head -n 5"


def test_run_cmd_handles_paths_with_spaces_and_redirection(tmp_path):
    src = tmp_path / "input dir" / "my file.txt"
    src.parent.mkdir()
    src.write_text("hello\n")
    dst = tmp_path / "output; rm -rf nothing.txt"
    util.run_cmd_via_subprocess(["cat", str(src), ">", str(dst)], check_files=[str(dst)], verbose=False)
    assert dst.read_text() == "hello\n"


def test_run_cmd_chained_commands(tmp_path):
    first, second = tmp_path / "a b.txt", tmp_path / "c d.txt"
    util.run_cmd_via_subprocess(
        ["touch", str(first), ";", "touch", str(second)], check_files=[str(first), str(second)], verbose=False
    )


def test_run_cmd_stdout_stderr_files(tmp_path):
    out, err = tmp_path / "out.txt", tmp_path / "err.txt"
    util.run_cmd_via_subprocess(
        ["bash", "-c", "echo to-stdout; echo to-stderr >&2"], stdout=str(out), stderr=str(err), verbose=False
    )
    assert out.read_text() == "to-stdout\n"
    assert err.read_text() == "to-stderr\n"


def test_run_cmd_relative_stdout_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    util.run_cmd_via_subprocess(["echo", "hi"], stdout="out.txt", stderr="err.txt", verbose=False)
    assert (tmp_path / "out.txt").read_text() == "hi\n"


def test_run_cmd_failure_raises_when_not_exiting():
    with pytest.raises(RuntimeError, match="Had an issue running"):
        util.run_cmd_via_subprocess(["bash", "-c", "exit 3"], verbose=False, exit_on_error=False)


def test_run_cmd_failure_exits_by_default():
    with pytest.raises(SystemExit):
        util.run_cmd_via_subprocess(["bash", "-c", "exit 3"], verbose=False)


def test_run_cmd_missing_check_file(tmp_path):
    with pytest.raises(RuntimeError, match="does not exist after running"):
        util.run_cmd_via_subprocess(
            ["true"], check_files=[str(tmp_path / "missing.txt")], verbose=False, exit_on_error=False
        )


def test_run_cmd_invalid_stdout_directory(tmp_path):
    with pytest.raises(RuntimeError):
        util.run_cmd_via_subprocess(
            ["true"], stdout=str(tmp_path / "missing_dir" / "out.txt"), stderr=str(tmp_path / "err.txt"),
            verbose=False, exit_on_error=False,
        )


def test_multi_process_safe(tmp_path):
    target = tmp_path / "with space.txt"
    assert util.multi_process_safe(["touch", str(target)]) == ("success", None)
    assert target.is_file()
    status, msg = util.multi_process_safe(["bash", "-c", "echo boom >&2; exit 1"])
    assert status == "error" and "boom" in msg


def test_run_cmd_simple(tmp_path):
    assert util.run_cmd_simple(["touch", str(tmp_path / "x y")]) is True
    assert util.run_cmd_simple(["false"]) is False


def test_touch_concat_and_gunzip_helpers(tmp_path):
    a, b, target = tmp_path / "a.txt", tmp_path / "b.txt", tmp_path / "t.txt"
    a.write_text("A\n")
    b.write_text("B\n")
    util.concat_files([str(a), str(b)], str(target))
    assert target.read_text() == "A\nB\n"
    util.concat_files([str(a)], str(target), append=True)
    assert target.read_text() == "A\nB\nA\n"

    touched = tmp_path / "check point.txt"
    util.touch_file(str(touched))
    util.touch_file(str(touched))
    assert touched.is_file() and touched.read_text() == ""

    gz = tmp_path / "seq.fasta.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(">s\nACGT\n")
    out = tmp_path / "seq.fasta"
    util.gunzip_file(str(gz), str(out))
    assert out.read_text() == ">s\nACGT\n"
    assert gz.is_file()


def _fake_diamond(tmp_path, recluster_writes_output):
    """ A stand-in `diamond` executable: cluster writes a 2-column table; recluster optionally does nothing. """
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    script = bin_dir / "diamond"
    recluster_body = (
        'printf "A\\tA\\nA\\tB\\nC\\tC\\n" > "$out"' if recluster_writes_output
        else 'echo "Recluster has been temporarily removed for v2.1.25. No action was taken." >&2'
    )
    script.write_text(
        "#!/bin/bash\n"
        "sub=$1; shift; out=''\n"
        'while [ $# -gt 0 ]; do if [ "$1" = "--out" ]; then out=$2; fi; shift; done\n'
        'if [ "$sub" = "cluster" ]; then printf "A\\tA\\nA\\tB\\nC\\tC\\n" > "$out"; fi\n'
        f'if [ "$sub" = "recluster" ]; then {recluster_body}; fi\n'
        "exit 0\n"
    )
    script.chmod(0o755)
    return bin_dir


@pytest.mark.parametrize("recluster_writes_output", [True, False])
def test_diamond_recluster_missing_output_falls_back(tmp_path, monkeypatch, capsys, recluster_writes_output):
    """ DIAMOND v2.2.x ships recluster as a no-op that exits 0 without writing --out (seen with v2.2.9). """
    bin_dir = _fake_diamond(tmp_path, recluster_writes_output)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    faa = tmp_path / "proteins.faa"
    faa.write_text(">A\nMKV\n>B\nMKV\n>C\nMKL\n")
    cluster_file = tmp_path / "clusters.tsv"
    refined_file = tmp_path / "refined_clusters.tsv"
    util.diamond_cluster_and_recluster(str(faa), str(cluster_file), str(refined_file), logging.getLogger("t"))
    assert refined_file.read_text() == cluster_file.read_text() == "A\tA\nA\tB\nC\tC\n"
    warned = "diamond recluster did not produce output" in capsys.readouterr().err
    assert warned is (not recluster_writes_output)
