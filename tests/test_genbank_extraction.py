"""
Tests for extracting gene-cluster regions from genome GenBanks (util.create_genbanks / util.create_genbank), which
fai uses to write one GenBank per gene cluster found in a target genome.
"""
import gzip
import random

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from zol import util


def _scaffold(name, length, rng, n_genes=30):
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    rec = SeqRecord(Seq(seq), id=name, name=name, description=f"{name} test scaffold",
                    annotations={"molecule_type": "DNA"})
    pos = 50
    for i in range(n_genes):
        glen = rng.randint(300, 1500)
        if pos + glen >= length:
            break
        strand = rng.choice([1, -1])
        rec.features.append(SeqFeature(SimpleLocation(pos, pos + glen, strand=strand), type="CDS",
                                       qualifiers={"locus_tag": [f"{name}_{i:04d}"], "translation": ["MK"]}))
        pos += glen + rng.randint(20, 400)
    # one spliced (joined) feature, to check compound locations are shifted correctly
    rec.features.append(SeqFeature(CompoundLocation([SimpleLocation(1000, 1300, strand=1),
                                                     SimpleLocation(1500, 1800, strand=1)]), type="CDS",
                                   qualifiers={"locus_tag": [f"{name}_join"], "translation": ["MK"]}))
    return rec


@pytest.fixture
def genome_gbk(tmp_path):
    """ A gzipped three-scaffold genome GenBank, like those in a prepTG database. """
    rng = random.Random(5)
    records = [_scaffold("contig_1", 40000, rng), _scaffold("contig_2", 25000, rng), _scaffold("contig_3", 15000, rng)]
    path = tmp_path / "genome.gbk.gz"
    with gzip.open(path, "wt") as handle:
        SeqIO.write(records, handle, "genbank")
    return str(path), {r.id: r for r in records}


def _regions(tmp_path, prefix):
    return [
        (str(tmp_path / f"{prefix}_a.gbk"), "contig_1", 900, 9000),
        (str(tmp_path / f"{prefix}_b.gbk"), "contig_1", 5000, 20000),   # overlaps region a on the same scaffold
        (str(tmp_path / f"{prefix}_c.gbk"), "contig_2", -300, 6000),    # starts before the scaffold
        (str(tmp_path / f"{prefix}_d.gbk"), "contig_1", 30000, 99999),  # runs past the scaffold end
        (str(tmp_path / f"{prefix}_e.gbk"), "contig_3", 2000, 7000),
    ]


def test_batch_matches_one_region_at_a_time(tmp_path, genome_gbk):
    path, _ = genome_gbk
    batch = _regions(tmp_path, "batch")
    single = _regions(tmp_path, "single")
    util.create_genbanks(path, batch)
    for region in single:
        util.create_genbank(path, *region)
    for (b, *_), (s, *_) in zip(batch, single):
        with open(b) as bh, open(s) as sh:
            assert bh.read() == sh.read()


def test_extracted_region_contents(tmp_path, genome_gbk):
    path, records = genome_gbk
    out = str(tmp_path / "region.gbk")
    start, end = 900, 9000
    util.create_genbanks(path, [(out, "contig_1", start, end)])
    rec = SeqIO.read(out, "genbank")
    original = records["contig_1"]
    assert str(rec.seq) == str(original.seq[start - 1:end])
    assert "Extracted region: contig_1:900..9000" in rec.annotations["comment"]

    expected = {}
    for f in original.features:
        if all(start <= p.start + 1 and p.end <= end for p in f.location.parts):
            expected[f.qualifiers["locus_tag"][0]] = [(p.start - start + 1, p.end - start + 1, p.strand)
                                                      for p in f.location.parts]
    got = {f.qualifiers["locus_tag"][0]: [(int(p.start), int(p.end), p.strand) for p in f.location.parts]
           for f in rec.features}
    assert got == expected
    assert "contig_1_join" in got and len(got["contig_1_join"]) == 2
    assert all(f.qualifiers["near_scaffold_edge"] in (["True"], ["False"]) for f in rec.features)


def test_overlapping_regions_do_not_affect_each_other(tmp_path, genome_gbk):
    """ The scaffold record is reused for every region on it, so extracting one region must not modify it. """
    path, _ = genome_gbk
    first = str(tmp_path / "first.gbk")
    overlap = str(tmp_path / "overlap.gbk")
    alone = str(tmp_path / "alone.gbk")
    util.create_genbanks(path, [(first, "contig_1", 900, 9000), (overlap, "contig_1", 5000, 20000)])
    util.create_genbanks(path, [(alone, "contig_1", 5000, 20000)])
    with open(overlap) as a, open(alone) as b:
        assert a.read() == b.read()


def test_genome_is_read_once_and_reading_stops_early(tmp_path, genome_gbk, monkeypatch):
    path, _ = genome_gbk
    parsed = []
    real_parse = SeqIO.parse

    def counting_parse(handle, fmt):
        parsed.append(None)
        for rec in real_parse(handle, fmt):
            parsed.append(rec.id)
            yield rec

    monkeypatch.setattr(util.SeqIO, "parse", counting_parse)
    regions = [(str(tmp_path / f"r{i}.gbk"), "contig_2", 100 + 2000 * i, 3000 + 2000 * i) for i in range(5)]
    regions.append((str(tmp_path / "r_first.gbk"), "contig_1", 100, 4000))
    util.create_genbanks(path, regions)
    assert parsed.count(None) == 1               # the file was opened and parsed once for all six regions
    assert [p for p in parsed if p] == ["contig_1", "contig_2"]  # stopped before contig_3, which was not needed
    assert all((tmp_path / f"r{i}.gbk").is_file() for i in range(5))


def test_missing_scaffold_fails_clearly(tmp_path, genome_gbk, capsys):
    path, _ = genome_gbk
    with pytest.raises(SystemExit):
        util.create_genbanks(path, [(str(tmp_path / "x.gbk"), "contig_404", 1, 100)])
    assert "Scaffold(s) contig_404 not found" in capsys.readouterr().err


def test_no_regions_is_a_no_op(tmp_path, genome_gbk):
    path, _ = genome_gbk
    util.create_genbanks(path, [])
    assert list(tmp_path.glob("*.gbk")) == []


def test_uncompressed_genbank(tmp_path, genome_gbk):
    path, records = genome_gbk
    plain = tmp_path / "genome.gbk"
    with gzip.open(path, "rt") as src, open(plain, "w") as dst:
        dst.write(src.read())
    out = str(tmp_path / "plain_region.gbk")
    util.create_genbanks(str(plain), [(out, "contig_3", 2000, 7000)])
    assert str(SeqIO.read(out, "genbank").seq) == str(records["contig_3"].seq[1999:7000])
