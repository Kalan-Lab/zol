"""
Static (PDF / PNG) figures for the zol suite, drawn with matplotlib, plus a numpy implementation of neighbor-joining.

These replace the R (ggplot2, gggenes, cowplot, ggtree, ape, ...) scripts zol previously wrote and ran via Rscript,
so R is no longer a dependency. Functions here raise exceptions on failure - callers handle reporting / exiting.

matplotlib is imported lazily (and with the non-interactive Agg backend) so that importing this module is cheap and
works on headless systems.
"""

import csv
import math
import textwrap
from collections import defaultdict

import numpy as np

NA_GREY = "#7f7f7f"  # ggplot2's default na.value ("grey50")
MM_TO_PT = 72.27 / 25.4  # ggplot2 text "size" is in mm


def _plt():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    return plt


def _to_float(value):
    try:
        f = float(str(value).strip().strip('"'))
    except (TypeError, ValueError):
        return math.nan
    return f if math.isfinite(f) else math.nan


def _read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ----------------------------------------------------------------------------------------------------------------------
# colour scales
# ----------------------------------------------------------------------------------------------------------------------
def _cmap(colors, name="zol"):
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name, list(colors))


def gradient2_colors(values, low, mid, high, midpoint=0.0, limits=None, na_color=NA_GREY):
    """
    Map values to colours like ggplot2's scale_fill_gradient2(): the midpoint gets the mid colour and the scale extends
    symmetrically to the value furthest from the midpoint. Values that are NaN, or outside `limits` if given, get
    na_color. Returns (list of hex colours, (vmin, vmax) of the shown range, colormap, normalizer).
    """
    from matplotlib.colors import Normalize, to_hex

    vals = np.asarray([_to_float(v) for v in values], dtype=float)
    if limits is not None:
        lo, hi = limits
        extent = max(abs(lo - midpoint), abs(hi - midpoint))
        shown = (lo, hi)
    else:
        finite = vals[np.isfinite(vals)]
        extent = float(np.max(np.abs(finite - midpoint))) if finite.size else 1.0
        shown = (float(finite.min()), float(finite.max())) if finite.size else (midpoint - 1.0, midpoint + 1.0)
    if extent == 0:
        extent = 1.0
    if shown[0] == shown[1]:
        shown = (shown[0] - extent, shown[1] + extent) if limits is None else shown
    cmap = _cmap([low, mid, high])
    norm = Normalize(vmin=midpoint - extent, vmax=midpoint + extent)
    colors = []
    for v in vals:
        if not np.isfinite(v) or (limits is not None and not (limits[0] <= v <= limits[1])):
            colors.append(na_color)
        else:
            colors.append(to_hex(cmap(norm(v))))
    return colors, shown, cmap, norm


def _colorbar(fig, cax, cmap, norm, shown, orientation="horizontal", label=None, label_size=10, ticks=None,
              ticklabels=None):
    from matplotlib.cm import ScalarMappable

    sm = ScalarMappable(norm=norm, cmap=cmap)
    lo, hi = shown
    if lo == hi:
        hi = lo + 1e-9
    cb = fig.colorbar(sm, cax=cax, orientation=orientation, boundaries=np.linspace(lo, hi, 257))
    cb.outline.set_linewidth(0.5)
    if ticks is None:
        from matplotlib.ticker import MaxNLocator

        ticks = [t for t in MaxNLocator(5).tick_values(lo, hi) if lo - 1e-12 <= t <= hi + 1e-12]
        cb.set_ticks(ticks)
    else:
        cb.set_ticks(ticks)
        if ticklabels is not None:
            cb.set_ticklabels(ticklabels)
    cb.ax.tick_params(labelsize=8, length=2)
    if label:
        cb.set_label(label, fontsize=label_size)
    return cb


# ----------------------------------------------------------------------------------------------------------------------
# gene arrows (replacement for gggenes::geom_gene_arrow)
# ----------------------------------------------------------------------------------------------------------------------
def draw_gene_arrows(ax, genes, y=0.0, body_height_mm=3.0, head_height_mm=4.0, head_length_mm=4.0,
                     edgecolor="black", linewidth=0.6):
    """
    Draw gene arrows on `ax`. genes: iterable of (start, end, forward, fill_colour). The arrow head is a fixed
    physical size (like gggenes); its length is capped at the gene length. The axes' x-limits and size must be final
    before calling this, since millimetres are converted to data units from them.
    """
    from matplotlib.patches import Polygon

    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x_per_mm = abs(x1 - x0) / max(bbox.width * 25.4, 1e-9)
    y_per_mm = abs(y1 - y0) / max(bbox.height * 25.4, 1e-9)
    bh = body_height_mm * y_per_mm / 2.0
    hh = head_height_mm * y_per_mm / 2.0
    for start, end, forward, color in genes:
        start, end = float(start), float(end)
        if end < start:
            start, end = end, start
        head = min(head_length_mm * x_per_mm, end - start)
        if forward:
            neck = end - head
            pts = [(start, y - bh), (neck, y - bh), (neck, y - hh), (end, y), (neck, y + hh), (neck, y + bh),
                   (start, y + bh)]
        else:
            neck = start + head
            pts = [(end, y - bh), (neck, y - bh), (neck, y - hh), (start, y), (neck, y + hh), (neck, y + bh),
                   (end, y + bh)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=edgecolor, linewidth=linewidth,
                             joinstyle="miter"))


def _clean_axis(ax, keep_y=False):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    if not keep_y:
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])


# ----------------------------------------------------------------------------------------------------------------------
# cgc
# ----------------------------------------------------------------------------------------------------------------------
def plot_cgc(rows, tracks, color_list, rel_heights, titles, out_file, height=7.0, width=7.0, squish_factor=0.0,
             bottom_spacing=0.1, show_copy_count_status=False, label_distance=0.025, label_text_size=2.0,
             legend_title_size=10.0, no_legend=False, png=False):
    """
    Collapsed gene cluster figure for cgc.

    rows: list of dicts (one per ortholog group in consensus order) with keys x_start, x_end, x_midpoint, direction
          (1 forward / 0 reverse), label, single_copy ("True"/"False") and one value per track name.
    tracks: track names, top to bottom; the last one colours the gene schematic, the others are bar tracks.
    color_list: [low, mid, high] colours per track. rel_heights: relative height per track (top to bottom).
    titles: display title per track.
    squish_factor: relative height of the gap between consecutive tracks (negative values make tracks overlap, as
                   with cowplot). bottom_spacing: relative height of extra space below the figure.
    """
    plt = _plt()
    fig = plt.figure(figsize=(width, height))
    xmin = min(float(r["x_start"]) for r in rows)
    xmax = max(float(r["x_end"]) for r in rows)
    pad = (xmax - xmin) * 0.02
    xlim = (xmin - pad, xmax + pad)

    show_sc = show_copy_count_status and any(r.get("single_copy") == "False" for r in rows)
    # cells from top to bottom: (kind, relative height)
    cells = []
    heights = [float(h) for h in rel_heights]
    for i, t in enumerate(tracks):
        if i > 0:
            cells.append(("gap", float(squish_factor)))
        if i < len(tracks) - 1:
            cells.append(("bar", heights[i], i))
        else:
            if show_sc:
                cells.append(("sc", 0.1 * heights[i], i))
                cells.append(("gap", float(squish_factor)))
                cells.append(("gene", 0.9 * heights[i], i))
            else:
                cells.append(("gene", heights[i], i))
    cells.append(("gap", float(bottom_spacing)))
    total = sum(c[1] for c in cells)
    if total <= 0:
        raise ValueError("The relative heights, squish factor and bottom spacing sum to a non-positive value.")

    left, right = 0.12, 0.98
    title_in = 0.25  # space reserved for track titles (inches)
    top_cursor = 1.0
    gene_ax = None
    gene_track = None
    for cell in cells:
        kind, rh = cell[0], cell[1]
        h = rh / total
        bottom = top_cursor - h
        if kind == "bar":
            i = cell[2]
            t = tracks[i]
            th = min(title_in / height, h * 0.4)
            ax = fig.add_axes([left, bottom, right - left, max(h - th, 1e-4)])
            vals = [_to_float(r.get(t)) for r in rows]
            cols, _, _, _ = gradient2_colors(vals, *color_list[i])
            for r, v, c in zip(rows, vals, cols):
                if not math.isfinite(v):
                    continue
                x0, x1 = float(r["x_start"]), float(r["x_end"])
                ax.bar(x0, v, width=x1 - x0, align="edge", color=c, edgecolor="black", linewidth=0.5)
            ax.set_xlim(*xlim)
            finite = [v for v in vals if math.isfinite(v)]
            lo = min([0.0] + finite)
            hi = max([0.0] + finite)
            if lo == hi:
                hi = lo + 1.0
            ax.set_ylim(lo - 0.05 * (hi - lo) if lo < 0 else lo, hi + 0.05 * (hi - lo))
            _clean_axis(ax, keep_y=True)
            ax.axhline(0, color="black", linewidth=0.5)
            ax.tick_params(axis="y", labelsize=8, length=2)
            fig.text(left, bottom + h - th * 0.1, titles[i], ha="left", va="top", fontsize=11)
        elif kind == "sc":
            ax = fig.add_axes([left, bottom, right - left, max(h, 1e-4)])
            xs = [float(r["x_midpoint"]) for r in rows if r.get("single_copy") == "False"]
            ax.scatter(xs, [0] * len(xs), color="#a33f37", s=8, zorder=3)
            ax.set_xlim(*xlim)
            ax.set_ylim(-1, 1)
            ax.axis("off")
        elif kind == "gene":
            gene_track = cell[2]
            legend_h = 0.0 if no_legend else min(0.7 / height, h * 0.45)
            gene_ax = fig.add_axes([left, bottom + legend_h, right - left, max(h - legend_h, 1e-4)])
            gene_ax.set_xlim(*xlim)
            gene_ax.set_ylim(0, 1)
            gene_ax.axis("off")
            gene_cell = (bottom, h, legend_h)
        top_cursor = bottom

    if gene_ax is not None:
        t = tracks[gene_track]
        cols, shown, cmap, norm = gradient2_colors([r.get(t) for r in rows], *color_list[gene_track])
        arrow_y = 0.72
        draw_gene_arrows(gene_ax, [(r["x_start"], r["x_end"], str(r["direction"]) == "1", c)
                                   for r, c in zip(rows, cols)], y=arrow_y)
        labels = [(float(r["x_midpoint"]), r.get("label", "")) for r in rows if r.get("label", "")]
        if labels:
            bbox = gene_ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            y_per_mm = 1.0 / max(bbox.height * 25.4, 1e-9)
            # label_distance is relative to the arrow position (default 0.025 ~= just below the arrow)
            label_y = arrow_y - 2.0 * y_per_mm - (label_distance / 0.025) * 1.0 * y_per_mm
            for x, lab in labels:
                gene_ax.text(x, label_y, lab, rotation=45, ha="right", va="top", rotation_mode="anchor",
                             fontsize=label_text_size * MM_TO_PT)
        if not no_legend:
            bottom, h, legend_h = gene_cell
            cax = fig.add_axes([0.5 - 0.2, bottom + legend_h * 0.35, 0.4, legend_h * 0.3])
            _colorbar(fig, cax, cmap, norm, shown, label=titles[gene_track], label_size=legend_title_size)
            cax.xaxis.set_label_position("top")

    fig.savefig(out_file, dpi=600 if png else 300)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------------
# generateSyntenicVisual.py
# ----------------------------------------------------------------------------------------------------------------------
def plot_syntenic_visual(track_file, out_file, height=3.0, width=7.0):
    """ Gene arrows coloured by an evolutionary statistic (limits -3 to 3, values outside shown in grey). """
    plt = _plt()
    rows = _read_tsv(track_file)
    if not rows:
        raise ValueError(f"No rows in {track_file}")
    cols, shown, cmap, norm = gradient2_colors([r["Metric"] for r in rows], "#e05c74", "#f2f2f2", "#2087b3",
                                               limits=(-3.0, 3.0))
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0.04, 0.35, 0.92, 0.55])
    xmin = min(float(r["Start"]) for r in rows)
    xmax = max(float(r["End"]) for r in rows)
    pad = (xmax - xmin) * 0.02
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(-1, 1)
    ax.axis("off")
    draw_gene_arrows(ax, [(r["Start"], r["End"], str(r["Direction"]) == "1", c) for r, c in zip(rows, cols)])
    for r in rows:
        if r.get("SC"):
            ax.text((float(r["Start"]) + float(r["End"])) / 2.0, 0, r["SC"], ha="center", va="center", fontsize=7)
    cax = fig.add_axes([0.3, 0.14, 0.4, 0.06])
    _colorbar(fig, cax, cmap, norm, shown, label="Metric", ticks=[-3, -2, -1, 0, 1, 2, 3],
              ticklabels=["", "-2", "", "0", "", "2", ""])
    cax.xaxis.set_label_position("top")
    fig.savefig(out_file)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------------
# fai
# ----------------------------------------------------------------------------------------------------------------------
def plot_table_pdf(tsv_file, pdf_file, rows_per_page=45, page_size=(10, 12)):
    """ Render a TSV table into a (multi-page) PDF. """
    plt = _plt()
    from matplotlib.backends.backend_pdf import PdfPages

    with open(tsv_file, newline="") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    header, body = table[0], table[1:]
    wrap = 70 if len(header) <= 2 else 40
    body = [["\n".join(textwrap.wrap(c, wrap)) or "" for c in row] for row in body]
    with PdfPages(pdf_file) as pdf:
        pages = [body[i:i + rows_per_page] for i in range(0, max(len(body), 1), rows_per_page)] or [[]]
        for page in pages:
            fig = plt.figure(figsize=page_size)
            ax = fig.add_axes([0.03, 0.03, 0.94, 0.94])
            ax.axis("off")
            tab = ax.table(cellText=page if page else [[""] * len(header)], colLabels=header, loc="upper center",
                           cellLoc="left", colLoc="left")
            tab.auto_set_font_size(False)
            tab.set_fontsize(8)
            tab.auto_set_column_width(list(range(len(header))))
            for (r, _), cell in tab.get_celld().items():
                cell.set_edgecolor("#d0d0d0")
                if r == 0:
                    cell.set_facecolor("#e8e8e8")
                    cell.set_text_props(weight="bold")
                else:
                    n_lines = max(1, cell.get_text().get_text().count("\n") + 1)
                    cell.set_height(cell.get_height() * n_lines)
            pdf.savefig(fig)
            plt.close(fig)


def plot_segments(seg_data_file, pdf_file, height=5.0, width=15.0):
    """
    One page per target genome: for each gene-cluster segment found in it, a bar per CDS showing the CDS to query
    length ratio, filled by percent identity (grey to red), outlined in black for key proteins.
    """
    plt = _plt()
    from matplotlib.backends.backend_pdf import PdfPages

    rows = _read_tsv(seg_data_file)
    by_sample = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_sample[r["sample"]][r["segment_title"]].append(r)
    cmap = _cmap(["#bebebe", "#ff0000"])
    from matplotlib.colors import Normalize

    norm = Normalize(0.0, 100.0)
    with PdfPages(pdf_file) as pdf:
        for sample in by_sample:  # order of first appearance, like R's unique()
            segments = by_sample[sample]
            titles = list(segments.keys())
            n_genes = [len({r["gene"] for r in segments[t]}) for t in titles]
            fig = plt.figure(figsize=(width, height))
            gs = fig.add_gridspec(2, len(titles), width_ratios=[max(n, 1) for n in n_genes],
                                  height_ratios=[1, 0.04], left=0.06, right=0.99, top=0.80, bottom=0.12,
                                  wspace=0.05, hspace=0.9)
            ymax = max([_to_float(r["sql_ratio"]) for r in rows if r["sample"] == sample] + [1.0])
            first_ax = None
            for j, t in enumerate(titles):
                ax = fig.add_subplot(gs[0, j], sharey=first_ax)
                first_ax = first_ax or ax
                genes = {}
                for r in segments[t]:
                    genes.setdefault(r["gene"], []).append(r)
                ordered = sorted(genes, key=lambda g: np.mean([_to_float(x["gene_order"]) for x in genes[g]]))
                for k, g in enumerate(ordered):
                    r = genes[g][-1]
                    iden = _to_float(r["identity"])
                    ax.bar(k, _to_float(r["sql_ratio"]), width=0.8,
                           color=cmap(norm(iden)) if math.isfinite(iden) else NA_GREY,
                           edgecolor="#000000" if r["key"] == "True" else "#ffffff", linewidth=1.5)
                ax.axhline(1.0, color="blue", linestyle=(0, (4, 3)), linewidth=0.8)
                ax.set_xticks(range(len(ordered)))
                ax.set_xticklabels(ordered, rotation=45, ha="right", fontsize=7)
                ax.set_xlim(-0.6, len(ordered) - 0.4)
                ax.set_ylim(0, ymax * 1.05)
                ax.grid(axis="y", color="#ebebeb", linewidth=0.6)
                ax.set_axisbelow(True)
                ax.set_title("\n".join(textwrap.wrap(t, 50)), fontsize=8, backgroundcolor="#d9d9d9")
                if j == 0:
                    ax.set_ylabel("CDS to Query Length Ratio")
                else:
                    ax.tick_params(axis="y", labelleft=False)
            fig.text(0.5, 0.01, "CDS Classifications", ha="center", va="bottom")
            cax = fig.add_axes([0.06, 0.935, 0.18, 0.025])
            _colorbar(fig, cax, cmap, norm, (0.0, 100.0), label="Percent identity", label_size=8, ticks=[0, 50, 100])
            cax.xaxis.set_label_position("top")
            fig.text(0.27, 0.947, "Outline: black = key protein, white = other", fontsize=8, va="center")
            pdf.savefig(fig)
            plt.close(fig)


def plot_tiny_aai(info_file, pdf_file):
    """ Scatter of AAI vs. proportion of query proteins found, coloured by mean syntenic correlation. """
    plt = _plt()
    rows = _read_tsv(info_file)
    x = np.array([_to_float(r["AAI"]) for r in rows], dtype=float)
    y = np.array([_to_float(r["Prop_Genes_Found"]) for r in rows], dtype=float)
    c = np.array([_to_float(r["Mean_Syntenic_Correlation"]) for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(10, 10))
    has_c = np.isfinite(c)
    if (~has_c).any():
        ax.scatter(x[~has_c], y[~has_c], color=NA_GREY, alpha=0.7, s=18, edgecolors="none", label="NA")
    if has_c.any():
        sc = ax.scatter(x[has_c], y[has_c], c=c[has_c], cmap=_cmap(["#e6ffbd", "#1754b0"]), alpha=0.7, s=18,
                        edgecolors="none")
        fig.colorbar(sc, ax=ax, shrink=0.6, label="Syntenic Correlation to Query")
    if (~has_c).any() and has_c.any():
        ax.legend(title="Syntenic Correlation to Query", loc="lower right", frameon=False)
    ax.set_xlabel("Average amino acid identity")
    ax.set_ylabel("Proportion of query proteins with match")
    ax.grid(color="#ebebeb", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.savefig(pdf_file)
    plt.close(fig)


def midpoint_root(tree):
    """
    Root an ete3 tree at the midpoint of its longest leaf-to-leaf path (as util.create_nj_tree does for prepTG's
    species trees), so user-provided unrooted / arbitrarily rooted trees are drawn consistently. Trees with fewer
    than 3 leaves, or already rooted at the midpoint, are returned unchanged.
    """
    if len(tree) < 3:
        return tree
    outgroup = tree.get_midpoint_outgroup()
    if outgroup is not None and outgroup is not tree:
        tree.set_outgroup(outgroup)
    # ete3's set_outgroup() puts the root halfway along the outgroup's branch - move it along that branch to the
    # true midpoint of the longest leaf-to-leaf path (which crosses the root)
    if len(tree.children) == 2:
        x_node, y_node = tree.children
        total = float(x_node.dist or 0.0) + float(y_node.dist or 0.0)

        def depth(node):
            return max((node.get_distance(leaf) for leaf in node), default=0.0) if not node.is_leaf() else 0.0

        dx, dy = depth(x_node), depth(y_node)
        x_len = min(max((dy + total - dx) / 2.0, 0.0), total)
        x_node.dist, y_node.dist = x_len, total - x_len
    return tree


def _tree_layout(tree):
    """
    Rectangular layout of an ete3 tree: returns (leaf names top to bottom, list of line segments [(x0,y0,x1,y1)]).
    x is the distance from the root; leaves are at y = 0, 1, 2, ...
    """
    pos = {}
    leaves = []
    segments = []
    root = tree
    xpos = {root: 0.0}
    for node in root.traverse("preorder"):
        if node is not root:
            xpos[node] = xpos[node.up] + max(float(node.dist or 0.0), 0.0)
    for node in root.traverse("postorder"):
        if node.is_leaf():
            pos[node] = float(len(leaves))
            leaves.append(node.name)
        else:
            ys = [pos[ch] for ch in node.children]
            pos[node] = (min(ys) + max(ys)) / 2.0
            segments.append((xpos[node], min(ys), xpos[node], max(ys)))
            for ch in node.children:
                segments.append((xpos[node], pos[ch], xpos[ch], pos[ch]))
    return leaves, segments


def plot_phylo_heatmap(tree_file, heatmap_data_file, pdf_file, height=7.0, width=10.0):
    """ Species tree (left) next to a heatmap of the best bitscore per query protein / OG in each genome. """
    from ete3 import Tree
    plt = _plt()
    from matplotlib.collections import LineCollection

    tree = midpoint_root(Tree(tree_file, format=1))
    leaves, segments = _tree_layout(tree)
    rows = _read_tsv(heatmap_data_file)
    ogs = sorted({r["query_prot_id"] for r in rows})
    og_index = {o: j for j, o in enumerate(ogs)}
    leaf_index = {name: i for i, name in enumerate(leaves)}
    mat = np.full((len(leaves), len(ogs)), np.nan)
    for r in rows:
        i = leaf_index.get(r["label"])
        v = _to_float(r["bitscore"])
        if i is None or not math.isfinite(v):
            continue
        j = og_index[r["query_prot_id"]]
        mat[i, j] = v if not math.isfinite(mat[i, j]) else max(mat[i, j], v)

    fig = plt.figure(figsize=(width, height))
    n = len(leaves)
    label_fs = max(2.0, min(8.0, 0.7 * 72.0 * (height * 0.8) / max(n, 1)))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.4, 1.0], left=0.03, right=0.80, top=0.97,
                          bottom=0.05 + min(0.3, 0.012 * max(len(o) for o in ogs) if ogs else 0.1), wspace=0.01)
    tax = fig.add_subplot(gs[0, 0])
    hax = fig.add_subplot(gs[0, 1], sharey=tax)
    tax.add_collection(LineCollection([[(a, b), (c, d)] for a, b, c, d in segments], colors="black",
                                      linewidths=0.6))
    xmax = max([s[2] for s in segments] + [1e-9])
    tax.set_xlim(-0.02 * xmax, xmax * 1.02)
    tax.set_ylim(n - 0.5, -0.5)
    tax.axis("off")

    cmap = _cmap(["#bebebe", "#000000"]).with_extremes(bad="#ffffff")
    hax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap=cmap, interpolation="nearest",
               extent=(-0.5, len(ogs) - 0.5, n - 0.5, -0.5))
    if len(ogs) <= 400:
        hax.set_xticks(range(len(ogs)))
        hax.set_xticklabels(ogs, rotation=90, fontsize=max(3.0, min(8.0, 500.0 / max(len(ogs), 1))))
    else:
        hax.set_xticks([])
    # white tile borders, like geom_tile(color="white"), when tiles are big enough to see them
    if n <= 300:
        hax.set_yticks(np.arange(-0.5, n, 1.0), minor=True)
        hax.set_xticks(np.arange(-0.5, len(ogs), 1.0), minor=True)
        hax.grid(which="minor", color="white", linewidth=0.5)
        hax.tick_params(which="minor", length=0)
    hax.yaxis.tick_right()
    if label_fs >= 2.5:
        hax.set_yticks(range(n))
        hax.set_yticklabels(leaves, fontsize=label_fs)
    else:
        hax.set_yticks([])
    hax.tick_params(axis="y", length=0, labelleft=False)
    for side in ("top", "right", "left"):
        hax.spines[side].set_visible(False)
    hax.set_xlabel("Query proteins / Homolog groups")
    fig.savefig(pdf_file)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------------
# salt / cgcg
# ----------------------------------------------------------------------------------------------------------------------
def plot_gc_vs_ribo_aai(input_data_file, pdf_file):
    """ Gene cluster AAI vs. ribosomal protein AAI, with a linear fit (and 95% CI) and the 1:1 line. """
    plt = _plt()
    rows = _read_tsv(input_data_file)
    x = np.array([_to_float(r["ribo_aai"]) for r in rows], dtype=float)
    y = np.array([_to_float(r["gc_aai"]) for r in rows], dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, alpha=0.7, s=12, color="black", edgecolors="none", zorder=3)
    if x.size:
        lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
        ax.plot([lo, hi], [lo, hi], linestyle=(0, (4, 3)), color="grey", linewidth=0.8, zorder=1)
    if x.size >= 3 and np.ptp(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 100)
        fit = slope * xs + intercept
        resid = y - (slope * x + intercept)
        dof = x.size - 2
        s = math.sqrt(float(np.sum(resid ** 2)) / dof) if dof > 0 else 0.0
        from scipy.stats import t as t_dist

        tval = float(t_dist.ppf(0.975, dof)) if dof > 0 else 0.0
        se = s * np.sqrt(1.0 / x.size + (xs - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
        ax.fill_between(xs, fit - tval * se, fit + tval * se, color="#999999", alpha=0.3, linewidth=0, zorder=2)
        ax.plot(xs, fit, color="red", linestyle=(0, (4, 3)), linewidth=1.0, zorder=4)
    ax.set_xlabel("Ribosomal Protein AAI")
    ax.set_ylabel("Gene Cluster AAI")
    ax.grid(color="#ebebeb", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(pdf_file)
    plt.close(fig)


def plot_color_legend(low_value, high_value, low_color, high_color, png_file):
    """ Stand-alone colour bar legend (PNG) for cgcg's node colouring. """
    plt = _plt()
    from matplotlib.colors import Normalize

    lo, hi = float(low_value), float(high_value)
    if lo == hi:
        hi = lo + 1e-9
    fig = plt.figure(figsize=(1.6, 4.0))
    cax = fig.add_axes([0.25, 0.08, 0.2, 0.84])
    _colorbar(fig, cax, _cmap([low_color, high_color]), Normalize(lo, hi), (lo, hi), orientation="vertical")
    cax.tick_params(labelsize=10)
    fig.savefig(png_file, dpi=600)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------------
# neighbor-joining (replacement for ape::nj)
# ----------------------------------------------------------------------------------------------------------------------
def neighbor_joining_newick(names, dist):
    """
    Build an unrooted neighbor-joining tree (Saitou & Nei 1987, as in ape::nj) from a symmetric distance matrix and
    return it as a Newick string. Runs in O(n^3) time with O(n^2) memory using vectorised numpy updates.
    """
    D = np.array(dist, dtype=float)
    n = len(names)
    if D.shape != (n, n):
        raise ValueError("Distance matrix shape does not match the number of names.")
    if not np.allclose(D, D.T, equal_nan=False):
        raise ValueError("Distance matrix is not symmetric.")
    if np.isnan(D).any():
        raise ValueError("Distance matrix contains missing values.")
    def bl(x):
        return repr(float(x))

    nodes = [str(x) for x in names]
    if n == 1:
        return f"{nodes[0]};"
    if n == 2:
        half = D[0, 1] / 2.0
        return f"({nodes[0]}:{bl(half)},{nodes[1]}:{bl(half)});"
    D = D.copy()
    np.fill_diagonal(D, 0.0)
    while len(nodes) > 3:
        m = len(nodes)
        r = D.sum(axis=1)
        Q = (m - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(np.argmin(Q)), m)
        if i > j:
            i, j = j, i
        dij = D[i, j]
        li = 0.5 * dij + (r[i] - r[j]) / (2.0 * (m - 2))
        lj = dij - li
        new_node = f"({nodes[i]}:{bl(li)},{nodes[j]}:{bl(lj)})"
        du = 0.5 * (D[i] + D[j] - dij)
        du[i] = 0.0
        D[i, :] = du
        D[:, i] = du
        D = np.delete(np.delete(D, j, axis=0), j, axis=1)
        nodes[i] = new_node
        del nodes[j]
    a = 0.5 * (D[0, 1] + D[0, 2] - D[1, 2])
    b = 0.5 * (D[0, 1] + D[1, 2] - D[0, 2])
    c = 0.5 * (D[0, 2] + D[1, 2] - D[0, 1])
    return f"({nodes[0]}:{bl(a)},{nodes[1]}:{bl(b)},{nodes[2]}:{bl(c)});"


def write_nj_tree_from_distance_file(input_dist_file, output_tree_file):
    """ Read the square, tab-delimited distance matrix written by util.create_nj_tree and write a NJ Newick tree. """
    with open(input_dist_file) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        names = header[1:]
        row_names = []
        mat = []
        for row in reader:
            if not row:
                continue
            row_names.append(row[0])
            mat.append([float(x) for x in row[1:]])
    if row_names != names:
        raise ValueError(f"Row and column names differ in distance matrix {input_dist_file}.")
    with open(output_tree_file, "w") as out:
        out.write(neighbor_joining_newick(names, mat) + "\n")


# ----------------------------------------------------------------------------------------------------------------------
# tree of gene clusters based on shared ortholog group content (for the interactive zol heatmap)
# ----------------------------------------------------------------------------------------------------------------------
MAX_UNIQUE_PROFILES_FOR_TREE = 12000  # a u x u float32 distance matrix: ~576 MB at this limit


def jaccard_distance_matrix(profiles, block=2048):
    """
    Square (u x u, float32) Jaccard distance matrix between boolean profiles (rows), computed in row blocks with
    matrix products so that no larger-than-output temporaries are created. Two empty profiles are at distance 0.
    """
    X = np.asarray(profiles, dtype=np.float32)
    u = X.shape[0]
    counts = X.sum(axis=1)
    D = np.empty((u, u), dtype=np.float32)
    for i in range(0, u, block):
        inter = X[i:i + block] @ X.T
        union = counts[i:i + block, None] + counts[None, :] - inter
        with np.errstate(divide="ignore", invalid="ignore"):
            blk = 1.0 - inter / union
        blk[union == 0] = 0.0
        D[i:i + block] = blk
    np.fill_diagonal(D, 0.0)
    return D


def average_linkage(D, sizes=None):
    """
    Average linkage (UPGMA) clustering of a square distance matrix, where item i initially stands for sizes[i]
    identical objects (default 1). This equals UPGMA over all the underlying objects, with identical objects merged
    first at height 0. Uses the nearest-neighbor chain algorithm with Lance-Williams updates: O(u^2) time, and it
    works in place on D (pass a copy to keep it).

    Returns a list of merges [a, b, height] in the order they were made (children always precede their parent),
    where ids < u are the items and id u + k is the node created by merge k.
    """
    D = np.asarray(D)
    u = D.shape[0]
    sizes = np.ones(u, dtype=np.float64) if sizes is None else np.asarray(sizes, dtype=np.float64).copy()
    if u < 2:
        return []
    np.fill_diagonal(D, np.inf)
    node_id = list(range(u))  # current tree node id held in each slot
    active = np.ones(u, dtype=bool)
    merges = []
    chain = []
    while len(merges) < u - 1:
        if not chain:
            chain.append(int(np.argmax(active)))
        a = chain[-1]
        row = D[a]
        b = int(np.argmin(row))
        if len(chain) > 1 and row[chain[-2]] <= row[b]:
            b = chain[-2]
        if len(chain) > 1 and b == chain[-2]:
            chain.pop()
            chain.pop()
            h = float(D[a, b])
            sa, sb = sizes[a], sizes[b]
            new_row = (sa * D[a] + sb * D[b]) / (sa + sb)
            new_row[a] = np.inf
            new_row[b] = np.inf
            new_row[~active] = np.inf
            # keep the merged cluster in slot a, retire slot b
            D[a, :] = new_row
            D[:, a] = new_row
            D[b, :] = np.inf
            D[:, b] = np.inf
            active[b] = False
            sizes[a] = sa + sb
            merges.append([node_id[a], node_id[b], h])
            node_id[a] = u + len(merges) - 1
        else:
            chain.append(b)
    # floating point rounding can make a parent minutely lower than a child - keep heights monotone
    heights = {}
    for k, m in enumerate(merges):
        m[2] = max([m[2]] + [heights[c] for c in m[:2] if c >= u])
        heights[u + k] = m[2]
    return merges


def og_content_tree(presence, max_unique=MAX_UNIQUE_PROFILES_FOR_TREE):
    """
    Average-linkage (UPGMA) tree of gene clusters based on the Jaccard distance between their sets of ortholog
    groups (OGs).

    presence: (n_clusters x n_ogs) boolean array.

    To scale to thousands of gene clusters, clusters with identical OG content are collapsed first (they are at
    distance 0 from each other) and the clustering runs on the u unique OG profiles, each weighted by how many gene
    clusters share it - which gives exactly the UPGMA tree of all n gene clusters. Memory is O(u^2) and time
    O(u^2 * n_ogs). Returns None if there are fewer than 2 gene clusters or more than `max_unique` unique profiles.

    Returns a dict with:
      "merges": list of [child_a, child_b, height] over all n gene clusters - ids < n are gene clusters and id n + k
                is the node created by merge k (children always precede their parent);
      "order": the leaf order (list of gene cluster indices, top to bottom);
      "unique_profiles": u.
    """
    P = np.asarray(presence, dtype=bool)
    n = P.shape[0]
    if n < 2:
        return None
    uniq, inverse = np.unique(P, axis=0, return_inverse=True)
    inverse = np.asarray(inverse).reshape(-1)
    u = uniq.shape[0]
    if u > max_unique:
        return None
    members = defaultdict(list)
    for c, g in enumerate(inverse):
        members[int(g)].append(c)

    merges = []
    next_id = n
    group_node = {}
    # gene clusters with identical profiles are joined at height 0
    for g in range(u):
        mem = members[g]
        node = mem[0]
        for c in mem[1:]:
            merges.append([node, c, 0.0])
            node = next_id
            next_id += 1
        group_node[g] = node

    if u >= 2:
        sizes = [len(members[g]) for g in range(u)]
        z_node = {g: group_node[g] for g in range(u)}
        for k, (a, b, h) in enumerate(average_linkage(jaccard_distance_matrix(uniq), sizes)):
            merges.append([z_node[a], z_node[b], h])
            z_node[u + k] = next_id
            next_id += 1

    # leaf order = left-to-right order of the tree
    kids = {n + k: (a, b) for k, (a, b, _) in enumerate(merges)}
    order, stack = [], [2 * n - 2]
    while stack:
        v = stack.pop()
        if v < n:
            order.append(v)
        else:
            stack.append(kids[v][1])
            stack.append(kids[v][0])
    return {"merges": merges, "order": order, "unique_profiles": int(u)}
