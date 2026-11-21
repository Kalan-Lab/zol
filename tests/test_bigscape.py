"""Tests for parsing BiG-SCAPE 1.x and 2.x clustering results and for zol-scape."""
import subprocess
import sys
from pathlib import Path

from zol import util

REPO_ROOT = Path(__file__).resolve().parents[1]

V1_CLUSTERING = (
    "#BGC Name\tFamily Number\n"
    "GCA_1.region001\t153\n"
    "GCA_2.region001\t153\n"
    "GCA_3.region002\t200\n"
)

V2_CLUSTERING = (
    "Record\tGBK\tRecord_Type\tRecord_Number\tCC\tFamily\n"
    "GCA_1.region001.gbk_region_1\tGCA_1.region001\tregion\t1\t12\tFAM_00034\n"
    "GCA_2.region001.gbk_region_1\tGCA_2.region001\tregion\t1\t12\tFAM_00034\n"
    "GCA_3.region002.gbk_region_2\tGCA_3.region002\tregion\t2\t13\tFAM_00040\n"
)

V2_PROTOCLUSTER_CLUSTERING = (
    "Record\tGBK\tRecord_Type\tRecord_Number\tCC\tFamily\n"
    "GCA_1.region001.gbk_protocluster_1\tGCA_1.region001\tprotocluster\t1\t1\tFAM_00001\n"
    "GCA_1.region001.gbk_protocluster_1_2\tGCA_1.region001\tprotocluster\t1_2\t1\tFAM_00001\n"
    "GCA_2.region001.gbk_protocluster_1\tGCA_2.region001\tprotocluster\t1\t1\tFAM_00001\n"
)

V2_BETA_CLUSTERING = (
    "GBK\tRecord_Type\tRecord_Number\tCC\tFamily\tFull_Name\n"
    "GCA_1.region001\tregion\t1\t12\tFAM_00034\tGCA_1.region001.gbk_region_1\n"
)


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)


def test_parse_v1(tmp_path):
    fmt, gcfs = util.parse_bigscape_clustering_file(_write(tmp_path / "NRPS_clustering_c0.30.tsv", V1_CLUSTERING))
    assert fmt == "v1"
    assert gcfs == {"153": ["GCA_1.region001", "GCA_2.region001"], "200": ["GCA_3.region002"]}


def test_parse_v2(tmp_path):
    fmt, gcfs = util.parse_bigscape_clustering_file(_write(tmp_path / "mix_clustering_c0.3.tsv", V2_CLUSTERING))
    assert fmt == "v2"
    assert gcfs == {"FAM_00034": ["GCA_1.region001", "GCA_2.region001"], "FAM_00040": ["GCA_3.region002"]}


def test_parse_v2_protocluster_dedupes_gbks(tmp_path):
    fmt, gcfs = util.parse_bigscape_clustering_file(
        _write(tmp_path / "mix_clustering_c0.3.tsv", V2_PROTOCLUSTER_CLUSTERING)
    )
    assert fmt == "v2"
    assert gcfs == {"FAM_00001": ["GCA_1.region001", "GCA_2.region001"]}


def test_parse_v2_beta(tmp_path):
    fmt, gcfs = util.parse_bigscape_clustering_file(_write(tmp_path / "mix_clustering_c0.3.tsv", V2_BETA_CLUSTERING))
    assert fmt == "v2"
    assert gcfs == {"FAM_00034": ["GCA_1.region001"]}


def test_parse_unknown_and_empty(tmp_path):
    assert util.parse_bigscape_clustering_file(_write(tmp_path / "a_clustering_x.tsv", "a\tb\tc\n1\t2\t3\n")) == \
        ("unknown", {})
    assert util.parse_bigscape_clustering_file(_write(tmp_path / "b_clustering_x.tsv", "")) == ("unknown", {})


def test_find_clustering_files(tmp_path):
    _write(tmp_path / "network_files" / "run" / "NRPS" / "NRPS_clustering_c0.30.tsv", V1_CLUSTERING)
    _write(tmp_path / "network_files" / "run" / "NRPS" / "NRPS_clans_0.30_0.70.tsv", "x")
    _write(tmp_path / "network_files" / "run" / "Network_Annotations_Full.tsv", "x")
    _write(tmp_path / "clustering.tsv", "gbk\trepresentative_gbk\n")
    found = util.find_bigscape_clustering_files(str(tmp_path))
    assert [Path(f).name for f in found] == ["NRPS_clustering_c0.30.tsv"]


def _make_bigscape_input(tmp_path):
    input_dir = tmp_path / "bigscape input"
    for bgc in ["GCA_1.region001", "GCA_2.region001", "GCA_3.region002"]:
        _write(input_dir / bgc.split(".")[0] / f"{bgc}.gbk", "LOCUS\n//\n")
    return input_dir


def _run_zol_scape_print_mode(tmp_path, results_dir, input_dir):
    outdir = tmp_path / "zol-scape out"
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "bin" / "zol-scape"), "-i", str(input_dir), "-r", str(results_dir),
         "-o", str(outdir), "-p", "-c", "2", "-z", "-s --domain-mode"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return (outdir / "zol.cmds").read_text().strip().splitlines()


def test_zol_scape_print_mode_v1(tmp_path):
    input_dir = _make_bigscape_input(tmp_path)
    results = tmp_path / "bigscape v1"
    _write(results / "network_files" / "2024-01-01_12-00-00_hybrids_glocal" / "NRPS" / "NRPS_clustering_c0.30.tsv",
           V1_CLUSTERING)
    cmds = _run_zol_scape_print_mode(tmp_path, results, input_dir)
    assert len(cmds) == 2
    assert all(cmd.startswith("zol -c 2 -s --domain-mode -o ") for cmd in cmds)
    assert any("GCA_1.region001.gbk" in cmd and "GCA_2.region001.gbk" in cmd for cmd in cmds)
    # paths containing spaces must be quoted
    assert all("'" in cmd for cmd in cmds)


def test_zol_scape_print_mode_v2(tmp_path):
    input_dir = _make_bigscape_input(tmp_path)
    results = tmp_path / "bigscape v2"
    _write(results / "output_files" / "2026-01-05_10-11-12_c0.3" / "mix" / "mix_clustering_c0.3.tsv", V2_CLUSTERING)
    _write(results / "output_files" / "2026-01-05_10-11-12_c0.3" / "record_annotations.tsv", "Record\tGBK\n")
    cmds = _run_zol_scape_print_mode(tmp_path, results, input_dir)
    assert len(cmds) == 2
    assert any("ClusterFile_0_FAM_00034" in cmd for cmd in cmds)
    assert any("GCA_3.region002.gbk" in cmd for cmd in cmds)
