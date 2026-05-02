"""
report.py
Generates ranked CSV output and interactive HTML report for AlleleSelect results.
HTML report includes:
    - Top 20 candidates table with sortable columns
    - Allele selectivity vs. accessibility scatter plot (Plotly)
    - mRNA secondary structure arc diagram with candidate window overlays
    - Summary statistics section
"""

import csv
import os
import json
from datetime import datetime


def save_csv(candidates: list, output_path: str) -> None:
    """
    Save ranked candidates to CSV with all scored fields.

    Columns: ASO_ID, ASO_sequence, length, mRNA_start_position, mRNA_end_position,
             mutation_position_in_ASO, delta_G_mutant, delta_G_wildtype,
             allele_selectivity_ratio, Tm_mutant, Tm_wildtype,
             mRNA_accessibility_score, off_target_count, splice_risk,
             recommended_gapmer_pattern, priority_rank
    """
    fieldnames = [
        "priority_rank",
        "ASO_ID",
        "ASO_sequence",
        "length",
        "mRNA_start_position",
        "mRNA_end_position",
        "mutation_position_in_ASO",
        "delta_G_mutant_kcal_mol",
        "delta_G_wildtype_kcal_mol",
        "allele_selectivity_ratio_kcal_mol",
        "Tm_mutant_C",
        "Tm_wildtype_C",
        "mRNA_accessibility_score",
        "off_target_count",
        "splice_risk",
        "ps_toxicity_flag",
        "recommended_gapmer_pattern",
        "top_candidate",
    ]

    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for rank, c in enumerate(candidates, 1):
            row = {
                "priority_rank": rank,
                "ASO_ID": c.get("ASO_ID", f"AS_{rank}"),
                "ASO_sequence": c.get("aso_seq", ""),
                "length": c.get("length", len(c.get("aso_seq", ""))),
                "mRNA_start_position": c.get("mRNA_start", ""),
                "mRNA_end_position": c.get("mRNA_end", ""),
                "mutation_position_in_ASO": c.get("mutation_pos_in_aso", ""),
                "delta_G_mutant_kcal_mol": round(c.get("dG_mutant", 0), 3),
                "delta_G_wildtype_kcal_mol": round(c.get("dG_wildtype", 0), 3),
                "allele_selectivity_ratio_kcal_mol": round(c.get("allele_selectivity_ratio", 0), 3),
                "Tm_mutant_C": round(c.get("Tm_mutant", 0), 1),
                "Tm_wildtype_C": round(c.get("Tm_wildtype", 0), 1),
                "mRNA_accessibility_score": round(c.get("accessibility_score", 0.5), 3),
                "off_target_count": c.get("off_target_count", -1),
                "splice_risk": c.get("splice_risk", "UNKNOWN"),
                "ps_toxicity_flag": c.get("ps_toxicity_flag", False),
                "recommended_gapmer_pattern": c.get("recommended_gapmer_pattern", ""),
                "top_candidate": c.get("top_candidate", False),
            }
            writer.writerow(row)

    print(f"CSV saved: {output_path}")


def save_html_report(
    candidates: list,
    variant_label: str,
    output_path: str,
    mrna_structure: dict = None,
) -> None:
    """
    Generate interactive HTML report with sortable table and Plotly scatter plot.

    Parameters
    ----------
    candidates : list of scored candidate dicts
    variant_label : str, e.g. "CACNA1A c.575G>A (R192Q)"
    output_path : str
    mrna_structure : dict with 'per_base_unpaired' and 'mfe_structure' (optional)
    """
    top_20 = candidates[:20]
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Prepare data for scatter plot
    scatter_data = []
    for rank, c in enumerate(candidates[:100], 1):
        asr = c.get("allele_selectivity_ratio", 0)
        acc = c.get("accessibility_score", 0.5)
        ot = c.get("off_target_count", -1)
        is_top = c.get("top_candidate", False)
        scatter_data.append({
            "rank": rank,
            "aso_id": c.get("ASO_ID", ""),
            "asr": round(asr, 3),
            "accessibility": round(acc, 3),
            "off_target": ot,
            "top": is_top,
            "length": c.get("length", 0),
            "splice_risk": c.get("splice_risk", "N"),
            "aso_seq": c.get("aso_seq", ""),
        })

    # Table rows HTML
    table_rows = ""
    for rank, c in enumerate(top_20, 1):
        is_top = c.get("top_candidate", False)
        row_class = 'class="top-candidate"' if is_top else ""
        asr = c.get("allele_selectivity_ratio", 0)
        acc = c.get("accessibility_score", 0.5)
        ot = c.get("off_target_count", -1)
        ot_str = str(ot) if ot >= 0 else "N/A"
        sr = c.get("splice_risk", "?")
        table_rows += f"""
        <tr {row_class}>
            <td>{rank}</td>
            <td style="font-family:monospace;font-size:0.85em">{c.get('aso_seq','')}</td>
            <td>{c.get('length','')}</td>
            <td>{c.get('mRNA_start','')}-{c.get('mRNA_end','')}</td>
            <td>{c.get('mutation_pos_in_aso','')}</td>
            <td>{round(c.get('dG_mutant',0),2)}</td>
            <td>{round(c.get('dG_wildtype',0),2)}</td>
            <td style="font-weight:bold;color:{'#0d7a7a' if asr < -1.5 else ('#444' if asr < -1.0 else '#c0392b')}">{round(asr,3)}</td>
            <td>{round(acc,3)}</td>
            <td>{ot_str}</td>
            <td style="color:{'#c0392b' if sr=='Y' else '#27ae60'}">{sr}</td>
            <td style="font-size:0.8em">{c.get('recommended_gapmer_pattern','')}</td>
        </tr>"""

    # Build full HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>AlleleSelect Report — {variant_label}</title>
<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
<style>
  body {{ font-family: 'Times New Roman', Times, serif; margin: 0; padding: 20px 40px;
         background: #fafafa; color: #1a1a1a; line-height: 1.5; }}
  h1 {{ font-size: 1.6em; border-bottom: 2px solid #0d7a7a; padding-bottom: 8px;
       color: #0d7a7a; }}
  h2 {{ font-size: 1.2em; color: #333; margin-top: 2em; }}
  .meta {{ color: #666; font-size: 0.9em; margin-bottom: 2em; }}
  .summary-grid {{ display: flex; gap: 20px; flex-wrap: wrap; margin: 1em 0 2em 0; }}
  .summary-card {{ background: white; border: 1px solid #ddd; border-radius: 6px;
                   padding: 14px 20px; min-width: 160px; text-align: center; }}
  .summary-card .val {{ font-size: 2em; font-weight: bold; color: #0d7a7a; }}
  .summary-card .lbl {{ font-size: 0.85em; color: #666; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85em; margin: 1em 0; }}
  th {{ background: #0d7a7a; color: white; padding: 8px 6px; text-align: left;
        cursor: pointer; user-select: none; }}
  td {{ padding: 6px 6px; border-bottom: 1px solid #eee; }}
  tr:hover td {{ background: #f0f9f9; }}
  tr.top-candidate td {{ background: #e6f7f7; }}
  tr.top-candidate:hover td {{ background: #d0f0f0; }}
  .plot-container {{ background: white; border: 1px solid #ddd; border-radius: 6px;
                     padding: 10px; margin: 1em 0 2em 0; }}
  .threshold-note {{ font-size: 0.85em; color: #555; margin: 0.5em 0; }}
  footer {{ margin-top: 3em; font-size: 0.8em; color: #888; border-top: 1px solid #ddd; padding-top: 1em; }}
</style>
</head>
<body>
<h1>AlleleSelect: {variant_label}</h1>
<div class="meta">
  Generated {now} | The Xiu Lab &mdash; <a href="https://thexiulab.org">thexiulab.org</a> |
  <a href="https://github.com/axshoe/alleleselect">github.com/axshoe/alleleselect</a>
</div>

<div class="summary-grid">
  <div class="summary-card"><div class="val">{len(candidates)}</div><div class="lbl">Total candidates</div></div>
  <div class="summary-card"><div class="val">{sum(1 for c in candidates if c.get('allele_selectivity_ratio',0) < -1.0)}</div><div class="lbl">ASR &lt; -1.0</div></div>
  <div class="summary-card"><div class="val">{sum(1 for c in candidates if c.get('top_candidate',False))}</div><div class="lbl">Top candidates</div></div>
  <div class="summary-card"><div class="val">{sum(1 for c in candidates[:50] if c.get('off_target_count',1) == 0)}</div><div class="lbl">Zero off-targets (top 50)</div></div>
</div>

<h2>Priority Space: Allele Selectivity vs. mRNA Accessibility</h2>
<p class="threshold-note">Optimal candidates (top-right, teal): high accessibility + strong allele selectivity + zero off-targets.</p>
<div class="plot-container" id="scatter-plot" style="height:480px"></div>

<h2>Top 20 Ranked Candidates</h2>
<p class="threshold-note">Teal rows: top candidates (ASR &lt; -1.5 kcal/mol, mutation at optimal position).</p>
<table id="candidates-table">
  <thead><tr>
    <th onclick="sortTable(0)">Rank</th>
    <th>ASO Sequence (5'→3')</th>
    <th onclick="sortTable(2)">Length</th>
    <th>mRNA Position</th>
    <th onclick="sortTable(4)">Mut Pos</th>
    <th onclick="sortTable(5)">ΔG mut</th>
    <th onclick="sortTable(6)">ΔG wt</th>
    <th onclick="sortTable(7)">ASR</th>
    <th onclick="sortTable(8)">Accessibility</th>
    <th onclick="sortTable(9)">Off-targets</th>
    <th>Splice Risk</th>
    <th>Gapmer Pattern</th>
  </tr></thead>
  <tbody>{table_rows}</tbody>
</table>

<script>
// Scatter plot
const scatter = {json.dumps(scatter_data)};
const teal = '#0d7a7a', orange = '#e67e22', red = '#c0392b', grey = '#999';

const colors = scatter.map(d =>
  d.off_target === 0 ? (d.top ? teal : '#1a8a8a') :
  d.off_target === -1 ? grey :
  d.off_target <= 2 ? orange : red
);

const sizes = scatter.map(d => d.top ? 14 : 8);

const trace = {{
  x: scatter.map(d => d.accessibility),
  y: scatter.map(d => d.asr),
  mode: 'markers',
  type: 'scatter',
  marker: {{ color: colors, size: sizes, opacity: 0.8,
             line: {{ width: scatter.map(d => d.top ? 2 : 0.5), color: '#333' }} }},
  text: scatter.map(d => `${{d.aso_id}}<br>Seq: ${{d.aso_seq}}<br>ASR: ${{d.asr}}<br>Access: ${{d.accessibility}}<br>Off-targets: ${{d.off_target}}`),
  hoverinfo: 'text',
}};

const layout = {{
  xaxis: {{ title: 'mRNA Accessibility Score (0–1)', range: [0, 1] }},
  yaxis: {{ title: 'Allele Selectivity Ratio (kcal/mol)' }},
  shapes: [
    {{ type: 'line', x0: 0, x1: 1, y0: -1.0, y1: -1.0, line: {{ dash: 'dash', color: '#999', width: 1 }} }},
    {{ type: 'line', x0: 0, x1: 1, y0: -1.5, y1: -1.5, line: {{ dash: 'dash', color: teal, width: 1.5 }} }},
    {{ type: 'line', x0: 0.65, x1: 0.65, y0: -6, y1: 2, line: {{ dash: 'dot', color: '#aaa', width: 1 }} }},
  ],
  annotations: [
    {{ x: 0.67, y: -1.55, text: 'Top candidate zone', showarrow: false,
       font: {{ size: 11, color: teal }}, xanchor: 'left' }},
    {{ x: 0.02, y: -1.05, text: 'ASR = -1.0 threshold', showarrow: false,
       font: {{ size: 10, color: '#888' }}, xanchor: 'left' }},
  ],
  margin: {{ t: 20, r: 20, b: 50, l: 60 }},
  paper_bgcolor: 'white', plot_bgcolor: '#fafafa',
  font: {{ family: 'Times New Roman' }},
}};

Plotly.newPlot('scatter-plot', [trace], layout, {{responsive: true}});

// Table sort
function sortTable(col) {{
  const table = document.getElementById('candidates-table');
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.rows);
  rows.sort((a, b) => {{
    const va = parseFloat(a.cells[col].innerText) || a.cells[col].innerText;
    const vb = parseFloat(b.cells[col].innerText) || b.cells[col].innerText;
    return va < vb ? -1 : va > vb ? 1 : 0;
  }});
  rows.forEach(r => tbody.appendChild(r));
}}
</script>

<footer>
AlleleSelect v1.0.0 | The Xiu Lab | <a href="https://thexiulab.org">thexiulab.org</a> |
<a href="https://github.com/axshoe/alleleselect">github.com/axshoe/alleleselect</a><br>
This tool is for research use only. All candidates require experimental validation.
Results should not be used for clinical decision-making.
</footer>
</body>
</html>"""

    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML report saved: {output_path}")
