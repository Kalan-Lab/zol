"""Tests for pyhmmer usage - these should pass for all supported pyhmmer versions (>= 0.10.10)."""
import logging
import os

import pytest

from zol import util, zol


def _read_table(path):
    with open(path) as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip() and not line.startswith("#")]


def test_to_str_handles_bytes_str_and_none():
    assert util.pyhmmer_to_str(b"abc") == "abc"
    assert util.pyhmmer_to_str("abc") == "abc"
    assert util.pyhmmer_to_str(None) == "NA"


class _FakeQuery:
    name = "QueryHMM"
    accession = None


class _FakeTopHitsNew:
    """Mimics TopHits of pyhmmer >= 0.11 (query object)."""
    query = _FakeQuery()


class _FakeTopHitsOld:
    """Mimics TopHits of pyhmmer 0.10.x (query_name / query_accession attributes, bytes)."""
    query_name = b"QueryHMM"
    query_accession = b"PF00001.1"


def test_query_name_and_accession_across_api_versions():
    assert util.pyhmmer_query_name(_FakeTopHitsNew()) == "QueryHMM"
    assert util.pyhmmer_query_accession(_FakeTopHitsNew()) == "NA"
    assert util.pyhmmer_query_name(_FakeTopHitsOld()) == "QueryHMM"
    assert util.pyhmmer_query_accession(_FakeTopHitsOld()) == "PF00001.1"


def test_search_pfam_domains_coordinates(hmm_test_data):
    hits = zol._search_pfam_domains(
        hmm_test_data["protein_faa"], hmm_test_data["hmm_file"], hmm_test_data["z"],
        "Domain Gathering 1e-5", cpus=1, return_coordinates=True,
    )
    assert set(hits) == {"p1", "p2"}
    p1 = sorted(hits["p1"], key=lambda x: x[1])
    assert [d[0] for d in p1] == ["FAM0", "FAM1"]
    assert (p1[0][1], p1[0][2]) == (1, 120)
    assert (p1[1][1], p1[1][2]) == (125, 244)
    assert [d[0] for d in hits["p2"]] == ["FAM2"]
    assert all(isinstance(d[0], str) for d in p1)


@pytest.mark.parametrize("pfam_params", ["Domain Gathering 1e-5", "Full None 1e-5"])
def test_run_pyhmmer_pfam(hmm_test_data, tmp_path, pfam_params):
    out = tmp_path / "pfam.txt"
    status = zol.run_pyhmmer([
        "pfam", hmm_test_data["hmm_file"], hmm_test_data["z"], hmm_test_data["protein_faa"], str(out), 1, pfam_params,
    ])
    assert status == ("success", None)
    rows = {(r[0], r[2]): r for r in _read_table(out)}
    assert set(rows) == {("FAM0", "p1"), ("FAM1", "p1"), ("FAM2", "p2")}
    # accession is reported, and coordinates of the (included) domains are reported
    assert rows[("FAM0", "p1")][1] == "FAM0.1"
    assert rows[("FAM0", "p1")][3:5] == ["1", "120"]
    assert rows[("FAM1", "p1")][3:5] == ["125", "244"]


def test_run_pyhmmer_non_pfam(hmm_test_data, tmp_path):
    out = tmp_path / "other.txt"
    status = zol.run_pyhmmer([
        "custom", hmm_test_data["hmm_file"], hmm_test_data["z"], hmm_test_data["protein_faa"], str(out), 1, None,
    ])
    assert status == ("success", None)
    rows = {(r[0], r[2]): r for r in _read_table(out)}
    assert set(rows) == {("FAM0", "p1"), ("FAM1", "p1"), ("FAM2", "p2")}
    assert rows[("FAM2", "p2")][3:5] == ["1", "120"]


def test_run_pyhmmer_empty_proteome(hmm_test_data, tmp_path):
    out = tmp_path / "empty_res.txt"
    status = zol.run_pyhmmer([
        "custom", hmm_test_data["hmm_file"], hmm_test_data["z"], str(hmm_test_data["dir"] / "empty.faa"), str(out),
        1, None,
    ])
    assert status == ("success", None)
    assert _read_table(out) == []


def test_vog_annotation_coverage(hmm_test_data, tmp_path):
    out = tmp_path / "vog.txt"
    util.run_pyhmmer_for_vo_gfor_salt(
        [hmm_test_data["hmm_file"], hmm_test_data["z"], hmm_test_data["protein_faa"], str(out), 1]
    )
    rows = _read_table(out)
    assert {(r[0], r[2]) for r in rows} == {("FAM0", "p1"), ("FAM1", "p1"), ("FAM2", "p2")}
    for row in rows:
        assert 0.9 <= float(row[6]) <= 1.0


def test_ribosomal_protein_detection(hmm_test_data, tmp_path, monkeypatch):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    data_dir = tmp_path / "db"
    data_dir.mkdir()
    (data_dir / "database_location_paths.txt").write_text(
        f"riboprots\tNA\t{hmm_test_data['hmm_file']}\tNA\n"
    )
    monkeypatch.setenv("ZOL_DATA_PATH", str(data_dir))

    rec = SeqRecord(Seq("A" * 900), id="contig1", annotations={"molecule_type": "DNA"})
    for prot in SeqIO.parse(hmm_test_data["protein_faa"], "fasta"):
        rec.features.append(SeqFeature(
            SimpleLocation(0, 30, strand=1), type="CDS",
            qualifiers={"locus_tag": [prot.id], "translation": [str(prot.seq)]},
        ))
    gbk = tmp_path / "genome.gbk"
    SeqIO.write(rec, str(gbk), "genbank")
    query_file = tmp_path / "queries.faa"
    query_file.write_text("")
    work_dir = str(tmp_path) + "/"

    util.run_pyhmmer_for_ribo_prots(str(gbk), str(query_file), work_dir, logging.getLogger("test"))

    found = {rec.id for rec in SeqIO.parse(str(query_file), "fasta")}
    assert found == {"Ribosomal_Protein|p1", "Ribosomal_Protein|p2"}
