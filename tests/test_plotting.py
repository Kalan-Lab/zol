"""
Tests for zol.plotting (matplotlib figures and neighbor-joining that replaced the former R scripts) and for the tree
shown next to the interactive zol heatmap.
"""
import json
import re
from pathlib import Path

import numpy as np
import pytest

from zol import heatmap_html, plotting

REPO_ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------------------------------------------------------
# neighbor-joining
# ----------------------------------------------------------------------------------------------------------------------
def _patristic(newick, names):
    from ete3 import Tree

    t = Tree(newick)
    leaves = {leaf.name: leaf for leaf in t}
    return np.array([[leaves[a].get_distance(leaves[b]) for b in names] for a in names])


def test_nj_recovers_additive_tree():
    # textbook example (Saitou & Nei / Wikipedia): additive, so NJ must reproduce the distances exactly
    names = ["a", "b", "c", "d", "e"]
    D = np.array([[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]], float)
    nwk = plotting.neighbor_joining_newick(names, D)
    assert "np.float" not in nwk
    assert np.allclose(_patristic(nwk, names), D)


def test_nj_random_additive_trees():
    from ete3 import Tree

    rng = np.random.default_rng(3)
    for n in (4, 9, 30):
        t = Tree()
        t.populate(n, names_library=[f"s{i}" for i in range(n)], random_branches=True)
        names = [leaf.name for leaf in t]
        D = _patristic(t.write(format=5), names)
        nwk = plotting.neighbor_joining_newick(names, D)
        got = Tree(nwk)
        assert got.robinson_foulds(t, unrooted_trees=True)[0] == 0
        assert np.allclose(_patristic(nwk, names), D, atol=1e-9)
        rng.shuffle(names)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_nj_tiny_inputs(n):
    names = [f"x{i}" for i in range(n)]
    D = np.ones((n, n)) - np.eye(n)
    nwk = plotting.neighbor_joining_newick(names, D)
    assert nwk.endswith(";") and all(x in nwk for x in names)
    if n > 1:
        assert np.allclose(_patristic(nwk, names), D)


def test_nj_rejects_bad_matrices():
    with pytest.raises(ValueError):
        plotting.neighbor_joining_newick(["a", "b"], [[0, 1], [2, 0]])
    with pytest.raises(ValueError):
        plotting.neighbor_joining_newick(["a", "b", "c"], [[0, 1], [1, 0]])


def test_nj_from_distance_file(tmp_path):
    names = ["g1", "g2", "g3", "g4"]
    D = np.array([[0, .1, .3, .4], [.1, 0, .3, .4], [.3, .3, 0, .2], [.4, .4, .2, 0]])
    dist_file = tmp_path / "d.tsv"
    with open(dist_file, "w") as fh:
        fh.write("sample\t" + "\t".join(names) + "\n")
        for i, n in enumerate(names):
            fh.write(n + "\t" + "\t".join(str(np.float64(x)) for x in D[i]) + "\n")
    out = tmp_path / "t.nwk"
    plotting.write_nj_tree_from_distance_file(str(dist_file), str(out))
    from ete3 import Tree

    t = Tree(str(out))
    assert sorted(leaf.name for leaf in t) == names
    # same midpoint rooting step util.create_nj_tree applies
    t.set_outgroup(t.get_midpoint_outgroup())


# ----------------------------------------------------------------------------------------------------------------------
# tree of gene clusters from shared ortholog group content
# ----------------------------------------------------------------------------------------------------------------------
def _check_tree(res, n):
    merges = res["merges"]
    assert len(merges) == n - 1
    seen = set()
    for k, (a, b, h) in enumerate(merges):
        for child in (a, b):
            assert child < n + k, "children must be created before their parent"
            assert child not in seen, "every node is merged exactly once"
            seen.add(child)
    assert sorted(res["order"]) == list(range(n))
    # the leaf order is the left-to-right order of the tree
    kids = {n + k: (a, b) for k, (a, b, _) in enumerate(merges)}
    order, stack = [], [2 * n - 2]
    while stack:
        v = stack.pop()
        if v < n:
            order.append(v)
        else:
            stack.extend([kids[v][1], kids[v][0]])
    assert order == res["order"]


def test_og_content_tree_structure_and_identical_profiles():
    presence = [
        [1, 1, 0, 0],
        [1, 1, 0, 0],   # identical to 0
        [1, 1, 1, 0],
        [0, 0, 1, 1],
        [1, 1, 0, 0],   # identical to 0
        [0, 0, 1, 1],   # identical to 3
    ]
    res = plotting.og_content_tree(presence)
    _check_tree(res, 6)
    assert res["unique_profiles"] == 3
    merges = res["merges"]
    # identical profiles are joined at height 0, before anything else
    zero = [m for m in merges if m[2] == 0.0]
    assert len(zero) == 3
    # UPGMA heights never decrease towards the root
    height = {6 + k: m[2] for k, m in enumerate(merges)}
    assert all(m[2] >= height.get(c, 0.0) for m in merges for c in m[:2])
    # the root height is the average Jaccard distance between the {0,1,2,4} and {3,5} groups
    # d({0,1,4} vs {3,5}) = 1 for 6 pairs, d({2} vs {3,5}) = 1 - 1/4 for 2 pairs
    assert merges[-1][2] == pytest.approx((6 * 1.0 + 2 * 0.75) / 8, rel=1e-5)
    # gene clusters with identical content are adjacent in the order
    pos = {c: i for i, c in enumerate(res["order"])}
    assert abs(pos[3] - pos[5]) == 1


def test_average_linkage_with_sizes_equals_full_upgma():
    """ Collapsing identical items into weighted ones must give exactly the UPGMA tree of all items. """
    from scipy.cluster.hierarchy import cophenet, linkage
    from scipy.spatial.distance import pdist, squareform

    rng = np.random.default_rng(11)
    U = rng.random((40, 5))
    counts = rng.integers(1, 6, 40)
    full = np.repeat(U, counts, axis=0)
    ref = squareform(cophenet(linkage(pdist(full), "average")))
    merges = plotting.average_linkage(squareform(pdist(U)), counts)
    members = {i: [i] for i in range(40)}
    C = np.zeros((40, 40))
    for k, (a, b, h) in enumerate(merges):
        for x in members[a]:
            for y in members[b]:
                C[x, y] = C[y, x] = h
        members[40 + k] = members.pop(a) + members.pop(b)
    grp = np.repeat(np.arange(40), counts)
    assert np.allclose(C[np.ix_(grp, grp)], ref)


def test_og_content_tree_scales_by_collapsing_identical_profiles():
    rng = np.random.default_rng(0)
    base = rng.random((50, 40)) < 0.5
    presence = base[rng.integers(0, 50, 4000)]
    res = plotting.og_content_tree(presence)
    assert res["unique_profiles"] <= 50
    _check_tree(res, 4000)


def test_og_content_tree_limits():
    assert plotting.og_content_tree([[1, 0]]) is None
    presence = np.eye(30, dtype=bool)
    assert plotting.og_content_tree(presence, max_unique=10) is None
    res = plotting.og_content_tree(np.ones((5, 3), dtype=bool))  # all identical
    _check_tree(res, 5)
    assert all(m[2] == 0.0 for m in res["merges"])


def test_heatmap_html_embeds_tree(tmp_path):
    presence = [[1, 1, 0], [1, 0, 1], [1, 1, 0], [0, 1, 1]]
    res = plotting.og_content_tree(presence)
    ogs = [{"id": f"OG{j}", "label": f"OG{j}", "info": [], "metrics": {"m": float(j)}} for j in range(3)]
    clusters = [{"id": f"gc{i}", "label": f"gc{i}"} for i in range(4)]
    cells = [[j, i, 1, [f"lt{i}{j}"]] for i in range(4) for j in range(3) if presence[i][j]]
    data = {"title": "t", "ogs": ogs, "clusters": clusters, "cells": cells, "metrics": ["m"], "default_metric": "m",
            "row_orders": {"Tree": res["order"], "Alphabetical": [0, 1, 2, 3]}, "default_order": "Tree",
            "tree": {"merges": res["merges"]}, "tree_order": "Tree"}
    out = tmp_path / "h.html"
    heatmap_html.write_heatmap_html(data, str(out))
    html = out.read_text()
    embedded = json.loads(re.search(r'<script type="application/json" id="zol-data">(.*?)</script>', html,
                                    re.S).group(1).replace("<\\/", "</"))
    assert embedded["tree"]["merges"] == res["merges"]
    assert 'id="showtree"' in html
    # also fine without a tree
    data.pop("tree")
    heatmap_html.write_heatmap_html(data, str(out))


# ----------------------------------------------------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------------------------------------------------
def _assert_pdf(path, min_pages=1):
    data = Path(path).read_bytes()
    assert data.startswith(b"%PDF")
    assert data.count(b"/Type /Page\n") + data.count(b"/Type /Page ") + data.count(b"/Type/Page") >= min_pages \
        or b"/Page" in data


def _write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return str(path)


def _cgc_rows():
    rows, start = [], 1.0
    for i, (length, cons, taj, sc) in enumerate([(900, 1.0, -1.2, "True"), (1500, 0.8, 0.4, "False"),
                                                  (300, 0.5, "NA", "True"), (1200, 0.2, 2.1, "True")]):
        end = start + length
        rows.append({"Ortholog Group (OG) ID": f"OG_{i}", "conservation": str(cons), "tajimas_d": str(taj),
                     "single_copy": sc, "x_start": str(start), "x_end": str(end),
                     "x_midpoint": str((start + end) / 2), "direction": str(i % 2), "y_start": "0",
                     "label": f"OG_{i}"})
        start = end + 100
    return rows


@pytest.mark.parametrize("png", [False, True])
def test_plot_cgc(tmp_path, png):
    out = tmp_path / ("cgc.png" if png else "cgc.pdf")
    plotting.plot_cgc(_cgc_rows(), ["tajimas_d", "conservation"], [["#bd3131", "#ffffff", "#548ce8"],
                                                                  ["#ffffff", "#878787", "#000000"]],
                      ["1", "1"], ["Tajima's D", "Conservation"], str(out), show_copy_count_status=True,
                      squish_factor=-0.1, png=png)
    head = out.read_bytes()[:8]
    assert head.startswith(b"\x89PNG") if png else head.startswith(b"%PDF")


def test_plot_cgc_single_track_no_legend(tmp_path):
    out = tmp_path / "cgc.pdf"
    plotting.plot_cgc(_cgc_rows(), ["conservation"], [["#ffffff", "#878787", "#000000"]], ["1"], ["Conservation"],
                      str(out), no_legend=True)
    _assert_pdf(out)


def test_gradient2_matches_ggplot_semantics():
    from matplotlib.colors import to_rgb

    def close(a, b):
        return np.allclose(to_rgb(a), to_rgb(b), atol=1 / 128)  # colormaps are quantised to 256 levels

    cols, shown, _, _ = plotting.gradient2_colors([0.0, 1.0, "NA", -0.5], "#ff0000", "#ffffff", "#0000ff")
    assert close(cols[0], "#ffffff")     # midpoint 0 -> mid colour
    assert close(cols[1], "#0000ff")     # furthest value from the midpoint -> high colour
    assert cols[2] == plotting.NA_GREY   # missing -> grey50
    assert shown == (-0.5, 1.0)
    cols, _, _, _ = plotting.gradient2_colors([5.0, 3.0], "#ff0000", "#ffffff", "#0000ff", limits=(-3, 3))
    assert cols[0] == plotting.NA_GREY and close(cols[1], "#0000ff")  # out of limits -> NA colour, like ggplot's default oob


def test_plot_syntenic_visual(tmp_path):
    track = _write_tsv(tmp_path / "Track.txt", ["OG", "Start", "End", "Direction", "SC", "Metric"],
                       [["OG_1", 1, 900, 1, "", -1.0], ["OG_2", 1100, 2000, 0, "//", 2.5], ["OG_3", 2200, 2600, 1,
                                                                                           "", "nan"]])
    out = tmp_path / "Plot.pdf"
    plotting.plot_syntenic_visual(track, str(out))
    _assert_pdf(out)


def test_fai_plots(tmp_path):
    seg = _write_tsv(tmp_path / "Bar_Plot_Input.txt",
                     ["sample", "segment_title", "gene", "key", "gene_order", "identity", "sql_ratio"],
                     [["g1", "g1 Segment_1 (Retained); 10 bp to Scaffold Start; 20 bp to Scaffold End", "A|1", "True",
                       100, 99.0, 1.0],
                      ["g1", "g1 Segment_1 (Retained); 10 bp to Scaffold Start; 20 bp to Scaffold End", "nov_aux|1",
                       "False", 50, 30.0, 0.5],
                      ["g1", "g1 Segment_2 (Filtered); 5 bp to Scaffold Start; 9 bp to Scaffold End", "B|1", "False",
                       900, 70.0, 1.3],
                      ["g2", "g2 Segment_1 (Retained); 1 bp to Scaffold Start; 2 bp to Scaffold End", "A|1", "True",
                       10, 100.0, 1.0]])
    names = _write_tsv(tmp_path / "Name_Mapping.txt", ["NR_ID", "CDS_Locus_Tags"],
                       [[f"P{i}", ", ".join(f"lt_{i}_{j}" for j in range(i % 4 + 1))] for i in range(60)])
    plotting.plot_segments(seg, str(tmp_path / "seg.pdf"), height=5, width=15)
    plotting.plot_table_pdf(names, str(tmp_path / "names.pdf"))
    _assert_pdf(tmp_path / "seg.pdf")
    _assert_pdf(tmp_path / "names.pdf")
    aai = _write_tsv(tmp_path / "aai.txt", ["Sample", "AAI", "Prop_Genes_Found", "Mean_Syntenic_Correlation"],
                     [["a", 71.0, 0.8, "nan"], ["b", 98.5, 0.9, 0.95], ["c", 90.0, 0.5, 0.3]])
    plotting.plot_tiny_aai(aai, str(tmp_path / "aai.pdf"))
    _assert_pdf(tmp_path / "aai.pdf")


def test_plot_phylo_heatmap(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text("((g1:0.1,g2:0.2):0.05,(g3:0.3,g4:0.1):0.02);\n")
    hm = _write_tsv(tmp_path / "hm.txt", ["label", "query_prot_id", "bitscore"],
                    [["g1", "OG1", 300.0], ["g1", "OG1", 350.0], ["g1", "OG2", 120.0], ["g2", "OG2", 80.0],
                     ["g3", "OG1", 10.0], ["g4", "OG1", "NA"], ["not_in_tree", "OG1", 5.0]])
    plotting.plot_phylo_heatmap(str(tree), hm, str(tmp_path / "phylo.pdf"), height=5, width=8)
    _assert_pdf(tmp_path / "phylo.pdf")
    leaves, segments = plotting._tree_layout(__import__("ete3").Tree(str(tree)))
    assert leaves == ["g1", "g2", "g3", "g4"] and len(segments) == 3 + 6


def test_phylo_heatmap_tree_is_midpoint_rooted():
    from ete3 import Tree

    # unrooted (trifurcating) tree whose longest path a..d (0.1 + 0.2 + 0.9 = 1.2) has its midpoint on d's branch
    t = plotting.midpoint_root(Tree("(a:0.1,b:0.1,(c:0.2,d:0.9):0.2);"))
    assert len(t.children) == 2
    root_to = {leaf.name: t.get_distance(leaf) for leaf in t}
    assert root_to["d"] == pytest.approx(0.6) and root_to["a"] == pytest.approx(0.6)
    # already midpoint-rooted / ultrametric trees are left as they are
    u = Tree("((a:1,b:1):1,(c:1.5,d:1.5):0.5);")
    before = u.write(format=1)
    assert plotting.midpoint_root(u).write(format=1) == before


def test_og_content_tree_is_rooted_at_its_midpoint():
    """ UPGMA trees are ultrametric, so their root already is the midpoint of the longest leaf-to-leaf path. """
    from ete3 import Tree

    rng = np.random.default_rng(4)
    presence = rng.random((60, 25)) < 0.4
    res = plotting.og_content_tree(presence)
    n = 60
    nodes = {i: Tree(name=str(i)) for i in range(n)}
    height = {i: 0.0 for i in range(n)}
    for k, (a, b, h) in enumerate(res["merges"]):
        parent = Tree()
        for c in (a, b):
            child = nodes.pop(c)
            child.dist = h - height[c]
            parent.add_child(child)
        nodes[n + k], height[n + k] = parent, h
    root = nodes[2 * n - 2]
    dists = [root.get_distance(leaf) for leaf in root]
    assert max(dists) == pytest.approx(min(dists))       # ultrametric
    assert max(dists) == pytest.approx(res["merges"][-1][2])
    # half of the longest leaf-to-leaf path is the root height, i.e. midpoint rooting would not move the root
    outgroup = root.get_midpoint_outgroup()
    assert root.get_distance(outgroup) == pytest.approx(0.0, abs=1e-9) or outgroup.up is root


def test_salt_and_cgcg_plots(tmp_path):
    rng = np.random.default_rng(1)
    rows = [[f"gc{i}", x, x * 0.9 + rng.normal(0, 2)] for i, x in enumerate(rng.uniform(60, 100, 20))]
    data = _write_tsv(tmp_path / "aai.txt", ["gene_cluster", "ribo_aai", "gc_aai"], rows)
    plotting.plot_gc_vs_ribo_aai(data, str(tmp_path / "salt.pdf"))
    _assert_pdf(tmp_path / "salt.pdf")
    plotting.plot_color_legend(0.0, 1.0, "#ffffff", "#000000", str(tmp_path / "legend.png"))
    assert (tmp_path / "legend.png").read_bytes().startswith(b"\x89PNG")


def test_no_r_left():
    """ zol no longer writes or runs R scripts anywhere. """
    hits = []
    for folder in ("bin", "scripts", "src"):
        for path in (REPO_ROOT / folder).rglob("*"):
            if path.is_file() and path.suffix not in (".pyc", ".so", ".o") and path.name != "plotting.py":
                try:
                    text = path.read_text()
                except UnicodeDecodeError:
                    continue
                if "Rscript" in text or "library(ggplot2)" in text:
                    hits.append(str(path.relative_to(REPO_ROOT)))
    assert hits == []
