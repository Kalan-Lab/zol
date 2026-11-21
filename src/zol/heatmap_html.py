"""
Interactive (D3-based) HTML version of the zol heatmap overview showing presence / copy-count of ortholog groups
across gene clusters.

The page is fully self-contained: data is embedded as JSON and D3 is inlined from the copy bundled with gravis (a
zol dependency) so that the HTML works offline (e.g. when opened on an HPC system). Only if that copy can not be
found is D3 loaded from a CDN instead.
"""

import importlib.util
import json
import os

D3_CDN_URL = "https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"


def load_bundled_d3_js():
    """
    Return the source of D3 v7 bundled with gravis, unwrapped from gravis' require/define wrapper so it simply
    defines a global `d3` object - or None if it can not be found / parsed.
    """
    try:
        spec = importlib.util.find_spec("gravis")
        if spec is None or spec.origin is None:
            return None
        d3_file = os.path.join(
            os.path.dirname(spec.origin), "_internal", "plotting", "third_party", "d3", "d3.v7.min.def.js"
        )
        if not os.path.isfile(d3_file):
            return None
        with open(d3_file, encoding="utf-8") as odf:
            text = odf.read()
        start = text.find('function(t){"use strict"')
        end = text.rfind("// ));")
        if start < 0 or end < 0 or end <= start:
            return None
        body = text[start:end].rstrip()
        return (
            "// D3 v7 - https://d3js.org - Copyright 2010-2021 Mike Bostock - BSD 3-Clause license\n"
            "var d3 = {};\n(" + body + ")(d3);\n"
        )
    except Exception:
        return None


def _script_safe(text):
    # prevent embedded content from closing the surrounding <script> element
    return text.replace("</", "<\\/")


def write_heatmap_html(data, html_file):
    """
    Write the interactive heatmap HTML page.

    data: dictionary with keys "title", "ogs", "clusters", "cells", "metrics", "default_metric", "row_orders",
          "default_order" and optionally "tree" / "tree_order" (see zol.plot_interactive_heatmap for how it is
          constructed). data["tree"]["merges"] is a scipy-linkage-like list of [child_a, child_b, height] over the
          gene clusters (ids < number of clusters are clusters; id n + k is the node made by merge k), drawn as a
          dendrogram left of the heatmap whenever the gene clusters are ordered by data["tree_order"].
    """
    d3_js = load_bundled_d3_js()
    if d3_js is not None:
        d3_tag = "<script>\n" + _script_safe(d3_js) + "\n</script>"
    else:
        d3_tag = '<script src="' + D3_CDN_URL + '"></script>'

    html = (
        _HTML_TEMPLATE.replace("@@TITLE@@", _escape_html(data.get("title", "zol heatmap")))
        .replace("@@D3@@", d3_tag)
        .replace("@@DATA@@", _script_safe(json.dumps(data)))
    )
    with open(html_file, "w", encoding="utf-8") as ohf:
        ohf.write(html)


def _escape_html(text):
    return (
        str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
    )


_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>@@TITLE@@</title>
<style>
  :root {
    --bg: #ffffff; --fg: #1f2933; --muted: #6b7280; --line: #e5e7eb;
    --absent: #eef0f3; --present: #889cbd; --multi: #3f5a8a; --bar: #374151; --accent: #d97706;
  }
  * { box-sizing: border-box; }
  body { margin: 0; padding: 16px 20px 40px; background: var(--bg); color: var(--fg);
         font: 13px/1.4 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
  h1 { font-size: 18px; margin: 0 0 2px; }
  .sub { color: var(--muted); margin-bottom: 12px; }
  .controls { display: flex; flex-wrap: wrap; gap: 10px 18px; align-items: center; padding: 10px 12px;
              border: 1px solid var(--line); border-radius: 8px; margin-bottom: 12px; background: #fafbfc; }
  .controls label { display: flex; align-items: center; gap: 6px; color: var(--muted); }
  .controls select, .controls input[type=search] { font: inherit; padding: 3px 6px; border: 1px solid #cbd2d9;
              border-radius: 5px; background: #fff; color: var(--fg); }
  .controls input[type=search] { width: 220px; }
  .controls input[type=range] { width: 90px; }
  .controls button { font: inherit; padding: 4px 10px; border: 1px solid #cbd2d9; border-radius: 5px;
              background: #fff; cursor: pointer; }
  .controls button:hover { background: #f0f2f5; }
  #matchinfo { color: var(--accent); min-width: 10ch; }
  #treectl input[type=range] { width: 70px; }
  .plot { display: flex; align-items: flex-start; }
  #left { flex: none; }
  #right { flex: 1 1 auto; overflow-x: auto; overflow-y: hidden; }
  svg { display: block; }
  .legend { display: flex; flex-wrap: wrap; gap: 6px 18px; margin-top: 10px; color: var(--muted); }
  .legend span.sw { display: inline-block; width: 14px; height: 14px; vertical-align: -3px; margin-right: 5px;
              border: 1px solid #d0d5dc; }
  #tooltip { position: fixed; pointer-events: none; z-index: 10; max-width: 460px; padding: 8px 10px;
             background: rgba(255,255,255,0.97); border: 1px solid #cbd2d9; border-radius: 6px;
             box-shadow: 0 4px 14px rgba(0,0,0,0.12); font-size: 12px; display: none; }
  #tooltip .t { font-weight: 600; margin-bottom: 4px; }
  #tooltip table { border-collapse: collapse; }
  #tooltip td { padding: 1px 0; vertical-align: top; }
  #tooltip td.k { color: var(--muted); padding-right: 10px; white-space: nowrap; }
  #tooltip .lts { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 11px; word-break: break-all; }
</style>
</head>
<body>
<h1>@@TITLE@@</h1>
<div class="sub" id="subtitle"></div>
<div class="controls">
  <label>Top track <select id="metric"></select></label>
  <label>Order gene clusters <select id="order"></select></label>
  <label id="treectl"><input type="checkbox" id="showtree" checked> Tree <input type="range" id="tw" min="40" max="400"></label>
  <label>Cell width <input type="range" id="cw" min="2" max="40"></label>
  <label>Cell height <input type="range" id="ch" min="2" max="40"></label>
  <label>Highlight <input type="search" id="search" placeholder="OG, gene cluster, or annotation"></label>
  <span id="matchinfo"></span>
  <button id="download" type="button">Download SVG</button>
</div>
<div class="plot"><div id="left"></div><div id="right"></div></div>
<div class="legend">
  <span><span class="sw" style="background:var(--absent)"></span>Absent</span>
  <span><span class="sw" style="background:var(--present)"></span>Present (single copy)</span>
  <span><span class="sw" style="background:var(--multi)"></span>Present in multiple copies (number shown when cells are large enough)</span>
  <span>Columns: ortholog groups (OGs) in consensus order &middot; hover for details</span>
  <span id="treelegend"></span>
</div>
<div id="tooltip"></div>

@@D3@@
<script type="application/json" id="zol-data">@@DATA@@</script>
<script>
(function () {
  "use strict";
  const DATA = JSON.parse(document.getElementById("zol-data").textContent);
  const ogs = DATA.ogs, clusters = DATA.clusters;
  const nO = ogs.length, nC = clusters.length;
  const COLORS = { absent: "#eef0f3", present: "#889cbd", multi: "#3f5a8a", bar: "#374151",
                   accent: "#d97706", fg: "#1f2933", muted: "#6b7280" };
  const FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif";
  const BAR_H = 110, GAP = 8;

  const cellMap = new Map();
  const clusterOgCount = new Array(nC).fill(0);
  DATA.cells.forEach(function (c) {
    cellMap.set(c[0] * nC + c[1], { copies: c[2], lts: c[3] });
    clusterOgCount[c[1]] += 1;
  });

  // ---------- tree ----------
  // node ids: 0..nC-1 are gene clusters (leaves), nC + k is the node created by merge k
  let treeKids = null, treeH = null, treeSize = null, treeMaxH = 0;
  if (DATA.tree && DATA.tree.merges && DATA.tree.merges.length === nC - 1) {
    const nN = 2 * nC - 1;
    treeKids = new Array(nN).fill(null); treeH = new Float64Array(nN); treeSize = new Int32Array(nN).fill(1);
    DATA.tree.merges.forEach(function (m, k) {
      const id = nC + k;
      treeKids[id] = [m[0], m[1]]; treeH[id] = m[2];
      treeSize[id] = treeSize[m[0]] + treeSize[m[1]];
      if (m[2] > treeMaxH) treeMaxH = m[2];
    });
    if (treeMaxH <= 0) treeMaxH = 1;
  }
  function cladeLeaves(id) {
    const out = [], stack = [id];
    while (stack.length) {
      const v = stack.pop();
      if (v < nC) out.push(v); else { stack.push(treeKids[v][1]); stack.push(treeKids[v][0]); }
    }
    return out;
  }

  d3.select("#subtitle").text(nC + " gene clusters × " + nO + " ortholog groups");

  // ---------- controls ----------
  const state = {
    metric: DATA.default_metric,
    order: DATA.default_order,
    cw: Math.max(6, Math.min(28, Math.floor((window.innerWidth - 320) / Math.max(nO, 1)))),
    ch: nC > 2000 ? 2 : nC > 600 ? 4 : nC > 150 ? 6 : nC > 60 ? 10 : nC > 25 ? 14 : 18,
    query: "",
    showTree: true,
    tw: Math.max(60, Math.min(220, Math.round(window.innerWidth * 0.12))),
    clade: null       // index of the selected tree node (n + k for merge k) or null
  };
  const metricSel = d3.select("#metric");
  metricSel.selectAll("option").data(DATA.metrics).join("option")
    .attr("value", function (d) { return d; }).text(function (d) { return d; });
  metricSel.property("value", state.metric)
    .on("change", function () { state.metric = this.value; render(true); });
  const orderSel = d3.select("#order");
  orderSel.selectAll("option").data(Object.keys(DATA.row_orders)).join("option")
    .attr("value", function (d) { return d; }).text(function (d) { return d; });
  orderSel.property("value", state.order)
    .on("change", function () { state.order = this.value; state.clade = null; render(true); });
  d3.select("#cw").property("value", state.cw)
    .on("input", function () { state.cw = +this.value; render(false); });
  d3.select("#ch").property("value", state.ch)
    .on("input", function () { state.ch = +this.value; render(false); });
  const TREE = DATA.tree && DATA.tree.merges && DATA.tree.merges.length === nC - 1 ? DATA.tree : null;
  if (!TREE) d3.select("#treectl").style("display", "none");
  d3.select("#showtree").on("change", function () { state.showTree = this.checked; render(false); });
  d3.select("#tw").property("value", state.tw)
    .on("input", function () { state.tw = +this.value; render(false); });
  if (TREE) {
    d3.select("#treelegend").text("Tree: " + (TREE.method || "average linkage") + " of " +
      (TREE.metric || "Jaccard distance of OG content") + " (shown when ordering by \u201c" + DATA.tree_order +
      "\u201d; click a branch point to highlight its gene clusters)");
  }
  d3.select("#search").on("input", function () { state.query = this.value.trim().toLowerCase(); applyHighlight(); });
  d3.select("#download").on("click", downloadSVG);

  // ---------- helpers ----------
  const ctx = document.createElement("canvas").getContext("2d");
  function textWidth(s, px) { ctx.font = px + "px " + FONT; return ctx.measureText(s).width; }
  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }
  const fmtSmall = d3.format(".3~g"), fmtBig = d3.format(",.0f");
  function fmt(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return "NA";
    return Math.abs(v) >= 1000 ? fmtBig(v) : fmtSmall(v);
  }
  function metricValue(o) {
    const v = ogs[o].metrics[state.metric];
    return (v === null || v === undefined) ? null : v;
  }

  // ---------- tooltip ----------
  const tip = d3.select("#tooltip");
  function showTip(html, ev) {
    tip.html(html).style("display", "block");
    const node = tip.node(), w = node.offsetWidth, h = node.offsetHeight;
    let x = ev.clientX + 14, y = ev.clientY + 14;
    if (x + w > window.innerWidth - 8) x = Math.max(8, ev.clientX - w - 14);
    if (y + h > window.innerHeight - 8) y = Math.max(8, ev.clientY - h - 14);
    tip.style("left", x + "px").style("top", y + "px");
  }
  function hideTip() { tip.style("display", "none"); }
  function rowsHtml(pairs) {
    return "<table>" + pairs.map(function (p) {
      return "<tr><td class='k'>" + esc(p[0]) + "</td><td>" + esc(p[1]) + "</td></tr>";
    }).join("") + "</table>";
  }
  function ogPairs(o) {
    const og = ogs[o];
    const pairs = og.info.slice();
    const shown = new Set(pairs.map(function (p) { return p[0]; }));
    if (!shown.has(state.metric)) pairs.push([state.metric, fmt(metricValue(o))]);
    return pairs;
  }
  function ogTip(o) { return "<div class='t'>" + esc(ogs[o].label) + "</div>" + rowsHtml(ogPairs(o)); }
  function cellTip(o, c) {
    const cell = cellMap.get(o * nC + c);
    let html = "<div class='t'>" + esc(ogs[o].label) + " &middot; " + esc(clusters[c].id) + "</div>";
    if (cell) {
      html += "<div>Present: <b>" + cell.copies + "</b> cop" + (cell.copies === 1 ? "y" : "ies") + "</div>";
      if (cell.lts.length) {
        const shownLts = cell.lts.slice(0, 12);
        html += "<div class='lts'>" + shownLts.map(esc).join("<br>") +
          (cell.lts.length > shownLts.length ? "<br>… +" + (cell.lts.length - shownLts.length) + " more" : "") +
          "</div>";
      }
    } else {
      html += "<div style='color:" + COLORS.muted + "'>Absent from this gene cluster</div>";
    }
    return html + "<hr style='border:none;border-top:1px solid #e5e7eb;margin:6px 0'>" + rowsHtml(ogPairs(o));
  }
  function clusterTip(c) {
    return "<div class='t'>" + esc(clusters[c].id) + "</div>" +
      rowsHtml([["OGs present", clusterOgCount[c] + " of " + nO]]);
  }

  // ---------- layout / drawing ----------
  const left = d3.select("#left").append("svg").attr("font-family", FONT);
  const right = d3.select("#right").append("svg").attr("font-family", FONT);
  const lAxis = left.append("g");
  const lMetric = left.append("text").attr("text-anchor", "middle").attr("font-size", 10).attr("fill", COLORS.muted);
  const lRows = left.append("g");
  const gTree = left.append("g");
  const treePath = gTree.append("path").attr("class", "tree").attr("fill", "none").attr("stroke", "#4b5563").attr("stroke-width", 1)
    .attr("shape-rendering", "crispEdges");
  const cladePath = gTree.append("path").attr("fill", "none").attr("stroke", COLORS.accent).attr("stroke-width", 2);
  const nodeMark = gTree.append("circle").attr("class", "hl").attr("r", 4).attr("fill", COLORS.accent)
    .attr("pointer-events", "none").style("display", "none");
  const gBars = right.append("g");
  const gZero = right.append("line").attr("stroke", "#9aa1ab");
  const gHeat = right.append("g");
  const heatBg = gHeat.append("rect").attr("fill", COLORS.absent);
  const gCells = gHeat.append("g");
  const gGrid = gHeat.append("g").attr("stroke", "#ffffff").attr("class", "grid");
  const gCounts = gHeat.append("g").attr("text-anchor", "middle").attr("fill", "#ffffff");
  const colHL = gHeat.append("rect").attr("class", "hl").attr("fill", "none").attr("stroke", COLORS.fg)
    .attr("stroke-width", 1).attr("pointer-events", "none").style("display", "none");
  const rowHL = gHeat.append("rect").attr("class", "hl").attr("fill", "none").attr("stroke", COLORS.fg)
    .attr("stroke-width", 1).attr("pointer-events", "none").style("display", "none");
  const gCols = right.append("g");
  const overlay = right.append("rect").attr("class", "overlay").attr("fill", "transparent");

  let L = {};       // current layout
  let rowPos = [];  // rowPos[clusterIndex] = row position
  let treeOn = false, nodeX = null, nodeY = null;
  function treeX(h) { return 4 + (L.treeW - 10) * (1 - h / treeMaxH); }
  function computeTreeCoords() {
    nodeX = new Float64Array(2 * nC - 1); nodeY = new Float64Array(2 * nC - 1);
    for (let c = 0; c < nC; c++) { nodeX[c] = treeX(0); nodeY[c] = L.heatTop + rowPos[c] * L.ch + L.ch / 2; }
    for (let k = 0; k < nC - 1; k++) {
      const id = nC + k, a = treeKids[id][0], b = treeKids[id][1];
      nodeX[id] = treeX(treeH[id]); nodeY[id] = (nodeY[a] + nodeY[b]) / 2;
    }
  }
  function treePathD(ids) {
    const parts = [];
    ids.forEach(function (id) {
      const a = treeKids[id][0], b = treeKids[id][1], x = nodeX[id];
      parts.push("M" + nodeX[a].toFixed(1) + "," + nodeY[a].toFixed(1) + "H" + x.toFixed(1) +
        "V" + nodeY[b].toFixed(1) + "H" + nodeX[b].toFixed(1));
    });
    return parts.join("");
  }
  function internalIds(root) {
    const out = [], stack = [root];
    while (stack.length) {
      const v = stack.pop();
      if (v >= nC) { out.push(v); stack.push(treeKids[v][0]); stack.push(treeKids[v][1]); }
    }
    return out;
  }
  function drawTree() {
    if (!treeOn) { treePath.attr("d", null); cladePath.attr("d", null); nodeMark.style("display", "none"); return; }
    computeTreeCoords();
    const all = [];
    for (let k = 0; k < nC - 1; k++) all.push(nC + k);
    treePath.attr("d", treePathD(all));
    cladePath.attr("d", state.clade !== null ? treePathD(internalIds(state.clade)) : null);
  }
  function nearestNode(px, py) {
    let best = -1, bestD = 64;  // within 8 px
    for (let k = 0; k < nC - 1; k++) {
      const id = nC + k, dx = nodeX[id] - px, dy = nodeY[id] - py, d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = id; }
    }
    return best;
  }

  function render(animate) {
    const cw = state.cw, ch = state.ch;
    const order = DATA.row_orders[state.order] || clusters.map(function (_, i) { return i; });
    rowPos = new Array(nC);
    order.forEach(function (c, i) { rowPos[c] = i; });

    const rowFont = Math.min(11, ch - 1), colFont = Math.min(11, cw - 1);
    const showRowLabels = rowFont >= 5, showColLabels = colFont >= 5;
    const rowLabW = showRowLabels ? d3.max(clusters, function (c) { return textWidth(c.label, rowFont); }) || 0 : 0;
    const colLabH = showColLabels ? d3.max(ogs, function (o) { return textWidth(o.label, colFont); }) || 0 : 0;
    treeOn = !!(treeKids && state.showTree && state.order === DATA.tree_order);
    const treeW = treeOn ? state.tw : 0;
    const leftW = Math.ceil(Math.max(rowLabW + 14, 64)) + treeW;
    const heatTop = BAR_H + GAP, heatW = nO * cw, heatH = nC * ch;
    const totalH = heatTop + heatH + 6 + colLabH + 8;
    L = { cw: cw, ch: ch, heatTop: heatTop, heatW: heatW, heatH: heatH, barTop: 8, treeW: treeW };

    left.attr("width", leftW).attr("height", totalH);
    right.attr("width", heatW + 12).attr("height", totalH);

    const t = animate ? d3.transition().duration(450) : null;
    function tr(sel) { return t ? sel.transition(t) : sel; }

    // top track
    const vals = ogs.map(function (_, o) { return metricValue(o); });
    const finite = vals.filter(function (v) { return v !== null; });
    let lo = Math.min(0, d3.min(finite) || 0), hi = Math.max(0, d3.max(finite) || 0);
    if (lo === hi) hi = lo + 1;
    const y = d3.scaleLinear().domain([lo, hi]).nice().range([BAR_H, L.barTop]);
    lAxis.attr("transform", "translate(" + (leftW - 4) + ",0)").call(d3.axisLeft(y).ticks(4).tickSizeOuter(0));
    lMetric.attr("transform", "translate(10," + ((BAR_H + L.barTop) / 2) + ") rotate(-90)")
      .text(state.metric.length > 22 ? state.metric.slice(0, 21) + "…" : state.metric);
    gZero.attr("x1", 0).attr("x2", heatW).attr("y1", y(0)).attr("y2", y(0));
    const bars = gBars.selectAll("rect").data(ogs.map(function (_, o) { return o; }), function (o) { return o; })
      .join("rect").attr("fill", COLORS.bar);
    bars.attr("x", function (o) { return o * cw + (cw > 3 ? 0.5 : 0); })
      .attr("width", Math.max(cw - (cw > 3 ? 1 : 0), 0.5));
    tr(bars).attr("y", function (o) { const v = vals[o]; return v === null ? y(0) : Math.min(y(v), y(0)); })
      .attr("height", function (o) { const v = vals[o]; return v === null ? 0 : Math.abs(y(v) - y(0)); });

    // heatmap
    gHeat.attr("transform", "translate(0," + heatTop + ")");
    heatBg.attr("width", heatW).attr("height", heatH);
    const cellData = DATA.cells.map(function (c) { return c; });
    const cells = gCells.selectAll("rect").data(cellData, function (c) { return c[0] + ":" + c[1]; })
      .join("rect").attr("fill", function (c) { return c[2] > 1 ? COLORS.multi : COLORS.present; });
    cells.attr("x", function (c) { return c[0] * cw; }).attr("width", cw).attr("height", ch);
    tr(cells).attr("y", function (c) { return rowPos[c[1]] * ch; });

    const vLines = cw >= 4 ? d3.range(1, nO) : [], hLines = ch >= 4 ? d3.range(1, nC) : [];
    gGrid.selectAll("line.v").data(vLines).join("line").attr("class", "v")
      .attr("x1", function (i) { return i * cw; }).attr("x2", function (i) { return i * cw; })
      .attr("y1", 0).attr("y2", heatH);
    gGrid.selectAll("line.h").data(hLines).join("line").attr("class", "h")
      .attr("y1", function (i) { return i * ch; }).attr("y2", function (i) { return i * ch; })
      .attr("x1", 0).attr("x2", heatW);

    const countData = (cw >= 12 && ch >= 10) ? cellData.filter(function (c) { return c[2] > 1; }) : [];
    const counts = gCounts.selectAll("text").data(countData, function (c) { return c[0] + ":" + c[1]; })
      .join("text").attr("font-size", Math.min(10, ch - 2)).attr("dominant-baseline", "central")
      .text(function (c) { return c[2]; });
    counts.attr("x", function (c) { return c[0] * cw + cw / 2; });
    tr(counts).attr("y", function (c) { return rowPos[c[1]] * ch + ch / 2; });

    // row labels (left) and column labels (bottom)
    const rl = lRows.selectAll("text").data(showRowLabels ? clusters.map(function (_, i) { return i; }) : [],
      function (c) { return c; })
      .join("text").attr("text-anchor", "end").attr("dominant-baseline", "central")
      .attr("x", leftW - 6).attr("font-size", rowFont).text(function (c) { return clusters[c].label; });
    tr(rl).attr("y", function (c) { return heatTop + rowPos[c] * ch + ch / 2; });
    drawTree();

    gCols.attr("transform", "translate(0," + (heatTop + heatH + 6) + ")");
    gCols.selectAll("text").data(showColLabels ? ogs.map(function (_, o) { return o; }) : [], function (o) { return o; })
      .join("text").attr("text-anchor", "end").attr("dominant-baseline", "central").attr("font-size", colFont)
      .attr("transform", function (o) { return "translate(" + (o * cw + cw / 2) + ",0) rotate(-90)"; })
      .text(function (o) { return ogs[o].label; });

    overlay.attr("width", heatW).attr("height", totalH);
    colHL.attr("width", cw).attr("height", heatH);
    rowHL.attr("width", heatW).attr("height", ch);
    applyHighlight();
  }

  // ---------- hover ----------
  function clearHover() {
    hideTip(); colHL.style("display", "none"); rowHL.style("display", "none"); nodeMark.style("display", "none");
    lRows.selectAll("text").attr("font-weight", null);
    gCols.selectAll("text").attr("font-weight", null);
  }
  overlay.on("mousemove", function (ev) {
    const p = d3.pointer(ev, right.node());
    const o = Math.floor(p[0] / L.cw);
    if (o < 0 || o >= nO) { clearHover(); return; }
    colHL.style("display", null).attr("x", o * L.cw).attr("y", 0);
    gCols.selectAll("text").attr("font-weight", function (d) { return d === o ? 700 : null; });
    const yIn = p[1] - L.heatTop;
    if (yIn >= 0 && yIn < L.heatH) {
      const r = Math.floor(yIn / L.ch);
      const c = rowPos.indexOf(r);
      rowHL.style("display", null).attr("x", 0).attr("y", r * L.ch);
      lRows.selectAll("text").attr("font-weight", function (d) { return d === c ? 700 : null; });
      showTip(cellTip(o, c), ev);
    } else {
      rowHL.style("display", "none");
      lRows.selectAll("text").attr("font-weight", null);
      showTip(ogTip(o), ev);
    }
  }).on("mouseleave", clearHover);
  function cladeTip(id) {
    const leaves = cladeLeaves(id);
    const ex = leaves.slice(0, 8).map(function (c) { return clusters[c].id; });
    return "<div class='t'>Clade of " + treeSize[id] + " gene clusters</div>" +
      rowsHtml([["Merge height", fmt(treeH[id]) + " (average Jaccard distance)"]]) +
      "<div class='lts'>" + ex.map(esc).join("<br>") +
      (leaves.length > ex.length ? "<br>… +" + (leaves.length - ex.length) + " more" : "") + "</div>" +
      "<div style='color:" + COLORS.muted + ";margin-top:4px'>Click to " +
      (state.clade === id ? "clear the highlight" : "highlight these gene clusters") + "</div>";
  }
  left.on("click", function (ev) {
    if (!treeOn) return;
    const p = d3.pointer(ev, left.node());
    if (p[0] > L.treeW) return;
    const id = nearestNode(p[0], p[1]);
    state.clade = (id < 0 || id === state.clade) ? null : id;
    drawTree(); applyHighlight();
    if (id >= 0) showTip(cladeTip(id), ev);
  });
  left.on("mousemove", function (ev) {
    const p = d3.pointer(ev, left.node());
    if (treeOn && p[0] <= L.treeW) {
      const id = nearestNode(p[0], p[1]);
      rowHL.style("display", "none"); lRows.selectAll("text").attr("font-weight", null);
      if (id < 0) { nodeMark.style("display", "none"); hideTip(); left.style("cursor", null); return; }
      nodeMark.style("display", null).attr("cx", nodeX[id]).attr("cy", nodeY[id]);
      left.style("cursor", "pointer");
      showTip(cladeTip(id), ev);
      return;
    }
    nodeMark.style("display", "none"); left.style("cursor", null);
    const yIn = p[1] - L.heatTop;
    if (yIn < 0 || yIn >= L.heatH) { clearHover(); return; }
    const r = Math.floor(yIn / L.ch), c = rowPos.indexOf(r);
    rowHL.style("display", null).attr("x", 0).attr("y", r * L.ch);
    lRows.selectAll("text").attr("font-weight", function (d) { return d === c ? 700 : null; });
    showTip(clusterTip(c), ev);
  }).on("mouseleave", clearHover);

  // ---------- search highlighting ----------
  function applyHighlight() {
    const q = state.query;
    let ogMatch = null, clMatch = null;
    if (q) {
      ogMatch = new Set(); clMatch = new Set();
      ogs.forEach(function (og, o) {
        if (og.label.toLowerCase().indexOf(q) >= 0 ||
            og.info.some(function (p) { return String(p[1]).toLowerCase().indexOf(q) >= 0; })) ogMatch.add(o);
      });
      clusters.forEach(function (c, i) { if (c.id.toLowerCase().indexOf(q) >= 0) clMatch.add(i); });
      d3.select("#matchinfo").text(ogMatch.size + " OGs, " + clMatch.size + " gene clusters match");
    } else {
      d3.select("#matchinfo").text("");
    }
    const dimOgs = ogMatch && ogMatch.size > 0;
    const cladeSet = (treeOn && state.clade !== null) ? new Set(cladeLeaves(state.clade)) : null;
    if (cladeSet) {
      clMatch = (clMatch && clMatch.size > 0)
        ? new Set(Array.from(clMatch).filter(function (c) { return cladeSet.has(c); })) : cladeSet;
      d3.select("#matchinfo").text((q ? d3.select("#matchinfo").text() + "; " : "") +
        cladeSet.size + " gene clusters in selected clade");
    }
    const dimCls = (clMatch && clMatch.size > 0) || !!cladeSet;
    function ogOp(o) { return dimOgs && !ogMatch.has(o) ? 0.2 : 1; }
    function clOp(c) { return dimCls && !clMatch.has(c) ? 0.25 : 1; }
    gBars.selectAll("rect").attr("opacity", ogOp)
      .attr("fill", function (o) { return dimOgs && ogMatch.has(o) ? COLORS.accent : COLORS.bar; });
    gCells.selectAll("rect").attr("opacity", function (c) { return Math.min(ogOp(c[0]), clOp(c[1])); });
    gCounts.selectAll("text").attr("opacity", function (c) { return Math.min(ogOp(c[0]), clOp(c[1])); });
    gCols.selectAll("text").attr("opacity", ogOp)
      .attr("fill", function (o) { return dimOgs && ogMatch.has(o) ? COLORS.accent : COLORS.fg; });
    lRows.selectAll("text").attr("opacity", clOp)
      .attr("fill", function (c) { return dimCls && clMatch.has(c) ? COLORS.accent : COLORS.fg; });
  }

  // ---------- export ----------
  function downloadSVG() {
    clearHover();
    const ns = "http://www.w3.org/2000/svg";
    const lNode = left.node().cloneNode(true), rNode = right.node().cloneNode(true);
    const lw = +lNode.getAttribute("width"), rw = +rNode.getAttribute("width");
    const h = Math.max(+lNode.getAttribute("height"), +rNode.getAttribute("height"));
    const out = document.createElementNS(ns, "svg");
    out.setAttribute("xmlns", ns);
    out.setAttribute("width", lw + rw); out.setAttribute("height", h);
    out.setAttribute("font-family", FONT);
    const bg = document.createElementNS(ns, "rect");
    bg.setAttribute("width", lw + rw); bg.setAttribute("height", h); bg.setAttribute("fill", "#ffffff");
    out.appendChild(bg);
    const gl = document.createElementNS(ns, "g"), gr = document.createElementNS(ns, "g");
    while (lNode.firstChild) gl.appendChild(lNode.firstChild);
    while (rNode.firstChild) gr.appendChild(rNode.firstChild);
    gr.setAttribute("transform", "translate(" + lw + ",0)");
    out.appendChild(gl); out.appendChild(gr);
    out.querySelectorAll(".overlay, .hl").forEach(function (n) { n.remove(); });
    const blob = new Blob([new XMLSerializer().serializeToString(out)], { type: "image/svg+xml" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = "zol_heatmap.svg";
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 1000);
  }

  render(false);
})();
</script>
</body>
</html>
"""
