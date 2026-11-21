"""
Tests for parallel execution and error reporting:
- util.robust_multiprocess_executor / util.robust_serial_executor (the worker contract used across zol)
- util.report_cause (short terminal message, full traceback in the log)
- fai's hmmlearn-based gene-cluster HMM, including use from a process pool.
"""
import logging
import multiprocessing
import os
import pickle
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from zol import fai, util

import parallel_workers as pw

TESTS_DIR = Path(__file__).resolve().parent
SRC_DIR = TESTS_DIR.parent / "src"


@pytest.fixture
def file_logger(tmp_path):
    """ A logger writing to a file, like the Progress.log zol creates (not the shared "task_logger"). """
    log_file = tmp_path / "Progress.log"
    logger = logging.getLogger(f"zol_test_{id(tmp_path)}")
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(log_file)
    logger.addHandler(handler)
    yield logger, log_file
    handler.close()
    logger.removeHandler(handler)


EXECUTORS = {
    "multiprocess": lambda worker, inputs, **kw: util.robust_multiprocess_executor(worker, inputs, 3, **kw),
    "serial": lambda worker, inputs, **kw: util.robust_serial_executor(worker, inputs, **kw),
}


# ----------------------------------------------------------------------------------------------------------------------
# Worker contract
# ----------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_successes_and_none_count_as_success(executor):
    inputs = [("ok", i) for i in range(10)] + [("none", None)] * 3
    summary = EXECUTORS[executor](pw.status_worker, inputs, description="ok items")
    assert summary == {"success_count": 13, "error_count": 0, "total_processed": 13}


@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_report_and_continue_keeps_terminal_output_short(executor, capsys, file_logger):
    logger, log_file = file_logger
    detail = "Sample S2 failed\n  Error type: KeyError\n  Traceback:\n    line 1\n    line 2"
    inputs = [("ok", 1), ("error", detail), ("ok", 3), ("error", "Sample S4 failed")]
    summary = EXECUTORS[executor](pw.status_worker, inputs, error_strategy="report_and_continue",
                                  log_object=logger, description="genome parsing")
    assert summary["success_count"] == 2 and summary["error_count"] == 2 and summary["total_processed"] == 4

    err = capsys.readouterr().err
    assert "Warning - error in genome parsing: Sample S2 failed" in err
    assert "Sample S4 failed" in err
    assert "Traceback" not in err  # multi-line detail only goes to the log
    assert "2 of 4 items failed during genome parsing; see the log for details." in err

    log_text = log_file.read_text()
    assert "Error type: KeyError" in log_text and "line 2" in log_text


@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_worker_exceptions_are_reported_not_fatal(executor, capsys, file_logger):
    """ A worker that raises used to abort the whole loop (dropping later results); now it is one error. """
    logger, log_file = file_logger
    inputs = [("ok", 1), ("raise", "bad GenBank record"), ("ok", 2), ("ok", 3)]
    summary = EXECUTORS[executor](pw.status_worker, inputs, log_object=logger, description="parsing")
    assert summary["success_count"] == 3 and summary["error_count"] == 1
    err = capsys.readouterr().err
    assert "ValueError: bad GenBank record" in err
    assert "Traceback (most recent call last)" in log_file.read_text()


@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_sys_exit_in_worker_is_reported(executor, capsys):
    inputs = [("ok", 1), ("exit", 1), ("exit", 0), ("ok", 2)]
    summary = EXECUTORS[executor](pw.status_worker, inputs, description="exiting")
    assert summary["success_count"] == 3 and summary["error_count"] == 1
    assert "Worker exited with status 1" in capsys.readouterr().err


@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_collect_errors(executor):
    inputs = [("error", "first"), ("ok", 1), ("raise", "second")]
    summary = EXECUTORS[executor](pw.status_worker, inputs, error_strategy="collect_errors")
    assert summary["error_count"] == 2
    assert len(summary["errors"]) == 2
    assert "first" in summary["errors"]
    assert any(e.startswith("ValueError: second") for e in summary["errors"])


@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_ignore_and_continue_is_silent(executor, capsys):
    summary = EXECUTORS[executor](pw.status_worker, [("error", "x"), ("ok", 1)], error_strategy="ignore_and_continue")
    assert summary["error_count"] == 1
    assert capsys.readouterr().err.count("Warning") == 0


@pytest.mark.parametrize("executor", sorted(EXECUTORS))
def test_report_and_stop_raises_worker_error(executor):
    inputs = [("ok", 1), ("error", "Tree building failed for OG_7\nlong detail"), ("ok", 2)]
    with pytest.raises(util.WorkerError) as excinfo:
        EXECUTORS[executor](pw.status_worker, inputs, error_strategy="report_and_stop", description="trees")
    assert "Tree building failed for OG_7" in str(excinfo.value)
    assert "long detail" not in str(excinfo.value)
    assert isinstance(excinfo.value, RuntimeError)  # callers catching RuntimeError keep working


def test_report_and_stop_does_not_wait_for_remaining_work():
    """ Stopping terminates the pool instead of waiting for slow, already-queued items. """
    inputs = [("error", "fail fast")] + [("sleep", 20)] * 6
    start = time.time()
    with pytest.raises(util.WorkerError):
        util.robust_multiprocess_executor(pw.status_worker, inputs, 2, error_strategy="report_and_stop")
    assert time.time() - start < 15


def test_pool_level_failure_counts_unreturned_items_as_errors(capsys):
    """ A result that cannot be sent back from the worker (unpicklable) is a pool-level failure. """
    inputs = [("unpicklable", None)] + [("ok", i) for i in range(4)]
    summary = util.robust_multiprocess_executor(pw.status_worker, inputs, 1, description="unpicklable")
    assert summary["total_processed"] == len(inputs)
    assert summary["error_count"] >= 1
    assert "critical error in unpicklable" in capsys.readouterr().err

    with pytest.raises(util.WorkerError):
        util.robust_multiprocess_executor(pw.status_worker, inputs, 1, error_strategy="report_and_stop")


def test_empty_inputs():
    assert util.robust_multiprocess_executor(pw.status_worker, [], 4)["total_processed"] == 0
    assert util.robust_serial_executor(pw.status_worker, [])["total_processed"] == 0


def test_work_is_spread_over_processes():
    summary = util.robust_multiprocess_executor(pw.pid_worker, list(range(20)), 2, error_strategy="collect_errors")
    assert summary["success_count"] == 20


@pytest.mark.parametrize("start_method", [m for m in ("fork", "spawn") if m in multiprocessing.get_all_start_methods()])
def test_sys_exit_in_worker_does_not_hang_pool(start_method, tmp_path):
    """
    Regression test: SystemExit inside a multiprocessing.Pool worker kills that process and its result never
    arrives, so a bare imap would wait forever. Run in a subprocess with a timeout so a hang fails the test.
    fai forces the "spawn" start method on macOS, so both start methods are checked.
    """
    script = textwrap.dedent(f"""
        import multiprocessing, sys
        sys.path[:0] = [{str(SRC_DIR)!r}, {str(TESTS_DIR)!r}]
        from zol import util
        import parallel_workers as pw
        if __name__ == "__main__":
            multiprocessing.set_start_method({start_method!r}, force=True)
            inputs = [("ok", 1), ("exit", 1), ("ok", 2), ("exit", 2), ("ok", 3)]
            summary = util.robust_multiprocess_executor(pw.status_worker, inputs, 2)
            print(summary["success_count"], summary["error_count"])
    """)
    script_path = tmp_path / "run_pool.py"
    script_path.write_text(script)
    proc = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == ["3", "2"]
    assert "Worker exited with status 1" in proc.stderr
    assert "worker message before exiting" in proc.stderr


def test_common_case_parallel_outputs_with_some_failures(tmp_path, capsys):
    """ Typical zol step: one output file per sample written in parallel, where a few samples fail. """
    out_dir = tmp_path / "outputs"
    out_dir.mkdir()
    inputs = [(str(out_dir / f"sample_{i}.txt"), f"sample {i}\n") for i in range(12)]
    inputs.insert(5, (str(tmp_path / "missing_dir" / "sample_bad.txt"), "x"))
    summary = util.robust_multiprocess_executor(pw.write_file_worker, inputs, 4, description="writing outputs")
    assert summary["success_count"] == 12 and summary["error_count"] == 1
    assert sorted(p.name for p in out_dir.iterdir()) == sorted(f"sample_{i}.txt" for i in range(12))
    assert (out_dir / "sample_7.txt").read_text() == "sample 7\n"
    err = capsys.readouterr().err
    assert err.count("Warning - error in writing outputs") == 1
    assert "sample_bad.txt" in err


def test_parallel_shell_commands_with_a_failure(tmp_path, capsys):
    """ multi_process_safe is the worker zol uses for running external programs in parallel. """
    cmds = [["echo", f"hello {i}", ">", str(tmp_path / f"out {i}.txt")] for i in range(6)]
    cmds.append(["ls", str(tmp_path / "does not exist")])
    summary = util.robust_multiprocess_executor(util.multi_process_safe, cmds, 3, description="running commands")
    assert summary["success_count"] == 6 and summary["error_count"] == 1
    for i in range(6):
        assert (tmp_path / f"out {i}.txt").read_text() == f"hello {i}\n"
    assert "Command failed" in capsys.readouterr().err


# ----------------------------------------------------------------------------------------------------------------------
# report_cause
# ----------------------------------------------------------------------------------------------------------------------

def test_report_cause_one_line_to_stderr_traceback_to_log(capsys, file_logger):
    logger, log_file = file_logger
    try:
        {}["missing_sample"]
    except Exception:
        sys.stderr.write("Issues with parsing the target genome listing.\n")
        util.report_cause(logger)
    err = capsys.readouterr().err
    assert err == "Issues with parsing the target genome listing.\n  Cause: KeyError: 'missing_sample'\n"
    log_text = log_file.read_text()
    assert "Cause: KeyError: 'missing_sample'" in log_text and "Traceback (most recent call last)" in log_text


def test_report_cause_without_log(capsys):
    task_logger = logging.getLogger("task_logger")
    saved = list(task_logger.handlers)
    for h in saved:
        task_logger.removeHandler(h)
    try:
        try:
            raise OSError("disk full")
        except OSError:
            util.report_cause()
        assert capsys.readouterr().err == "  Cause: OSError: disk full\n"
        try:
            raise OSError("disk full")
        except OSError:
            util.report_cause(fallback_traceback=True)
        err = capsys.readouterr().err
        assert "Cause: OSError: disk full" in err and "Traceback (most recent call last)" in err
    finally:
        for h in saved:
            task_logger.addHandler(h)


def test_report_cause_quiet_for_failed_validation_checks(capsys):
    try:
        raise AssertionError()
    except AssertionError:
        util.report_cause()
    util.report_cause()  # outside an except block: no-op
    assert capsys.readouterr().err == ""


def test_first_line():
    assert util.first_line("\n\n  first  \nsecond") == "first"
    assert util.first_line("") == ""
    assert len(util.first_line("x" * 1000)) == 300


# ----------------------------------------------------------------------------------------------------------------------
# fai gene-cluster HMM (hmmlearn)
# ----------------------------------------------------------------------------------------------------------------------

DEFAULTS = dict(gc_emission_prob_with_hit=0.95, gc_emission_prob_without_hit=0.2,
                gc_to_gc_transition_prob=0.9, bg_to_bg_transition_prob=0.9)


def default_model(hgs=("OG_1", "OG_2", "OG_3", "OG_4")):
    return fai.build_gene_cluster_hmm(list(hgs), **DEFAULTS)


def test_hmm_finds_cluster_in_background():
    bg = ["background"] * 10
    cluster = ["OG_1", "OG_2", "background", "OG_3", "OG_4"]
    states = fai.predict_gene_cluster_states(default_model(), bg + cluster + bg)
    assert states == [fai.BG_STATE] * 10 + [fai.GC_STATE] * 5 + [fai.BG_STATE] * 10


def test_hmm_edge_cases():
    model = default_model()
    assert fai.predict_gene_cluster_states(model, []) == []
    assert fai.predict_gene_cluster_states(model, ["OG_1"]) == [fai.GC_STATE]
    assert fai.predict_gene_cluster_states(model, ["background"]) == [fai.BG_STATE]
    assert model.labels == ["background", "OG_1", "OG_2", "OG_3", "OG_4"]
    with pytest.raises(ValueError, match="OG_99"):
        fai.predict_gene_cluster_states(model, ["OG_99"])
    with pytest.raises(ValueError, match="gc_to_gc_transition_prob"):
        fai.build_gene_cluster_hmm(["OG_1"], 0.95, 0.2, 1.5, 0.9)


def test_hmm_ortholog_group_identity_does_not_matter():
    """ All ortholog groups share one emission probability, so predictions depend only on hit / no hit. """
    model = default_model()
    a = ["background", "OG_1", "OG_1", "background", "OG_1", "background", "background"]
    b = ["background", "OG_2", "OG_4", "background", "OG_3", "background", "background"]
    assert fai.predict_gene_cluster_states(model, a) == fai.predict_gene_cluster_states(model, b)


def test_hmm_pickles_and_predicts_identically_in_pool():
    """ fai sends the model to every worker; results from a pool must match serial predictions. """
    import random

    model = default_model()
    restored = pickle.loads(pickle.dumps(model))
    rng = random.Random(1)
    scaffolds = [[rng.choice(model.labels) if rng.random() < 0.4 else "background" for _ in range(rng.randint(1, 150))]
                 for _ in range(40)]
    serial = [fai.predict_gene_cluster_states(model, s) for s in scaffolds]
    assert serial == [fai.predict_gene_cluster_states(restored, s) for s in scaffolds]
    with multiprocessing.Pool(3) as pool:
        parallel = pool.map(pw.predict_worker, [(model, s) for s in scaffolds])
    assert parallel == serial


def test_hmm_matches_pomegranate_if_installed():
    """ Optional check against the previous pomegranate model (only runs where pomegranate is installed). """
    pytest.importorskip("pomegranate")
    import numpy as np
    from pomegranate.distributions import Categorical
    from pomegranate.hmm import DenseHMM

    rng = np.random.default_rng(7)
    for ge, be, gt, bt in [(0.95, 0.2, 0.9, 0.9), (0.8, 0.3, 0.7, 0.95)]:
        hgs = [f"OG_{i}" for i in range(5)]
        labels = ["background"] + hgs
        gc_cat = Categorical([[be] + [ge] * len(hgs)])
        bg_cat = Categorical([[1 - be] + [1 - ge] * len(hgs)])
        old = DenseHMM()
        old.add_distributions([gc_cat, bg_cat])
        for a, b, p in [(old.start, gc_cat, 0.5), (old.start, bg_cat, 0.5), (gc_cat, old.end, 0.5),
                        (bg_cat, old.end, 0.5), (gc_cat, gc_cat, gt), (gc_cat, bg_cat, 1 - gt),
                        (bg_cat, gc_cat, 1 - bt), (bg_cat, bg_cat, bt)]:
            old.add_edge(a, b, p)
        new = fai.build_gene_cluster_hmm(hgs, ge, be, gt, bt)
        for _ in range(50):
            seq = [labels[int(rng.integers(1, 6))] if rng.random() < 0.3 else "background"
                   for _ in range(int(rng.integers(2, 120)))]
            expected = old.predict(np.array([[[labels.index(h)] for h in seq]]))[0].numpy().tolist()
            assert fai.predict_gene_cluster_states(new, seq) == expected
