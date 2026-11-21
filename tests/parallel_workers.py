"""
Worker functions for tests/test_parallel.py.

They live in their own module so multiprocessing can import them in child processes (macOS uses "spawn").
Each worker takes a single item, like the workers zol passes to util.robust_multiprocess_executor.
"""
import os
import sys
import time


def status_worker(item):
    """ item = (kind, value). Mirrors the ('success'|'error', ...) contract used by zol's workers. """
    kind, value = item
    if kind == "ok":
        return ("success", value)
    if kind == "none":
        return None
    if kind == "error":
        return ("error", value)
    if kind == "raise":
        raise ValueError(value)
    if kind == "exit":
        sys.stderr.write(f"worker message before exiting: {value}\n")
        sys.exit(value)
    if kind == "sleep":
        time.sleep(value)
        return ("success", None)
    if kind == "unpicklable":
        return ("success", lambda: None)
    raise AssertionError(f"unknown kind {kind}")


def write_file_worker(item):
    """ item = (path, text). Writes a per-item output file, failing for a path in a missing directory. """
    path, text = item
    try:
        with open(path, "w") as handle:
            handle.write(text)
        return ("success", None)
    except Exception as e:
        return ("error", f"Could not write {path}: {e}")


def predict_worker(item):
    """ item = (model, hgs_ordered). Predicts gene-cluster states, as fai's workers do with a shared model. """
    from zol import fai

    model, hgs_ordered = item
    return fai.predict_gene_cluster_states(model, hgs_ordered)


def pid_worker(_item):
    return ("success", os.getpid())
