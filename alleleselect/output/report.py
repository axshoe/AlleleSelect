"""
alleleselect/output/report.py  (v2)

Writes candidates.csv and report.html.
v2: adds composite_score, snp_pos_in_aso, snp_pos_score, snp_region,
    tox_summary, tox_serious, tox_warning, tox_flags to CSV output.
"""

from __future__ import annotations
import csv
import os
import html as html_module
from typing import List, Dict, Any


# ─── CSV column order ─────────────────────────────────────────────────────────
# v1 columns preserved in original order; v2 columns appended at end.
CSV_FIELDS = [
    # --- v1 columns ---
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
    # --- v2 columns ---
    "composite_score",
    "snp_pos_in_aso",
    "snp_pos_score",
    "snp_region",
    "tox_summary",
    "tox_serious",
    "tox_warning",
    "tox_flags",
]

# Map from internal candidate dict keys → CSV column names
# (only needed where they differ)
_KEY_MAP = {
    "aso_seq":                    "ASO_sequence",
    "ASO_ID":                     "ASO_ID",
    "mRNA_start":                 "mRNA_start_position",
    "mRNA_end":                   "mRNA_end_position",
    "mutation_pos_in_aso":        "mutation_position_in_ASO",
    "delta_G_mutant":             "delta_G_mutant_kcal_mol",
    "delta_G_wildtype":           "delta_G_wildtype_kcal_mol",
    "allele_selectivity_ratio":   "allele_selectivity_ratio_kcal_mol",
    "Tm_mutant":                  "Tm_mutant_C",
    "Tm_wildtype":                "Tm_wildtype_C",
    "accessibility_score":        "mRNA_accessibility_score",
    "off_target_count":           "off_target_count",
    "splice_risk":                "splice_risk",
    "ps_toxicity_flag":           "ps_toxicity_flag",
    "recommended_gapmer_pattern": "recommended_gapmer_pattern",
    "top_candidate":              "top_candidate",
    # v2
    "composite_score":            "composite_score",
    "snp_pos_in_aso":             "snp_pos_in_aso",
    "snp_pos_score":              "snp_pos_score",
    "snp_region":                 "snp_region",
    "tox_summary":                "tox_summary",
    "tox_serious":                "tox_serious",
    "tox_warning":                "tox_warning",
    "tox_flags":                  "tox_flags",
}


def _candidate_to_row(rank: int, c: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten a candidate dict into a flat CSV row dict."""
    row: Dict[str, Any] = {}

    # Map all known keys
    for src_key, csv_col in _KEY_MAP.items():
        if src_key in c:
            row[csv_col] = c[src_key]

    # Fields that are already named correctly in the candidate dict
    for field in CSV_FIELDS:
        if field in c and field not in row:
            row[field] = c[field]

    # Derived / positional
    row["priority_rank"] = rank
    row["length"] = len(c.get("aso_seq", "")) or c.get("length", "")

    # Ensure all CSV_FIELDS have a value (default empty string)
    for f in CSV_FIELDS:
        if f not in row:
            row[f] = ""

    return row


def save_csv(candidates: List[Dict[str, Any]], path: str) -> None:
    """Write all candidates to a CSV file with v1 + v2 columns."""
    os.makedirs(os.path.dirname(path) if os.path.dirname(path) else ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for rank, c in enumerate(candidates, 1):
            row = _candidate_to_row(rank, c)
            writer.writerow(row)
    print(f"CSV saved: {path}")


# ─── HTML report ──────────────────────────────────────────────────────────────

def _fmt(val: Any, decimals: int = 3) -> str:
    if val == "" or val is None:
        return ""
    try:
        f = float(val)
        return f"{f:.{decimals}f}"
    except (ValueError, TypeError):
        return str(val)


def save_html_report(
    candidates: List[Dict[str, Any]],
    variant_label: str,
    path: str,
) -> None:
    """Write an interactive HTML report with sortable table."""
    os.makedirs(os.path.dirname(path) if os.path.dirname(path) else ".", exist_ok=True)

    top5 = candidates[:5]

    def esc(s: Any) -> str:
        return html_module.escape(str(s))

    rows_html = []
    for rank, c in enumerate(candidates, 1):
        asr = c.get("allele_selectivity_ratio", c.get("allele_selectivity_ratio_kcal_mol", ""))
        acc = c.get("accessibility_score", "")
        ot  = c.get("off_target_count", "")
        comp = c.get("composite_score", "")
        snp_pos  = c.get("snp_pos_in_aso", "")
        snp_scr  = c.get("snp_pos_score", "")
        snp_reg  = c.get("snp_region", "")
        tox      = c.get("tox_summary", "")
        seq      = c.get("aso_seq", c.get("ASO_sequence", ""))
        aso_id   = c.get("ASO_ID", "")
        length   = c.get("length", len(seq) if seq else "")
        gapmer   = c.get("recommended_gapmer_pattern", "")
        splice   = c.get("splice_risk", "")

        tox_class = ""
        if tox and str(tox).startswith("FAIL"):
            tox_class = "tox-fail"
        elif tox and str(tox).startswith("WARN"):
            tox_class = "tox-warn"

        ot_str = str(ot) if ot not in ("", None, -1) else ("unscreened" if ot == -1 else "")
        ot_class = "ot-hits" if (isinstance(ot, int) and ot > 0) else ""

        top_style = ' class="top-row"' if rank <= 5 else ""

        rows_html.append(f"""
        <tr{top_style}>
            <td class="center">{rank}</td>
            <td><code>{esc(aso_id)}</code></td>
            <td><code class="seq">{esc(seq)}</code></td>
            <td class="center">{esc(length)}</td>
            <td class="center num">{_fmt(asr)}</td>
            <td class="center num">{_fmt(acc)}</td>
            <td class="center num">{_fmt(comp)}</td>
            <td class="center">{esc(snp_pos)}</td>
            <td class="center">{_fmt(snp_scr)}</td>
            <td class="center">{esc(snp_reg)}</td>
            <td class="center {ot_class}">{esc(ot_str)}</td>
            <td class="center">{esc(splice)}</td>
            <td class="{tox_class}">{esc(tox)}</td>
            <td class="small">{esc(gapmer)}</td>
        </tr>""")

    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AlleleSelect v2 — {esc(variant_label)}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          margin: 0; padding: 20px; background: #f8f9fa; color: #2c2c2c; }}
  h1   {{ color: #0d7a7a; font-size: 1.5rem; margin-bottom: 4px; }}
  .sub {{ color: #666; font-size: 0.9rem; margin-bottom: 20px; }}
  .badge {{ display: inline-block; background: #0d7a7a; color: white;
            border-radius: 4px; padding: 2px 8px; font-size: 0.75rem;
            margin-left: 8px; vertical-align: middle; }}
  table {{ width: 100%; border-collapse: collapse; background: white;
           box-shadow: 0 1px 4px rgba(0,0,0,.08); border-radius: 8px;
           overflow: hidden; font-size: 0.82rem; }}
  th    {{ background: #0d7a7a; color: white; padding: 8px 10px;
           text-align: left; cursor: pointer; user-select: none; white-space: nowrap; }}
  th:hover {{ background: #0b6868; }}
  td    {{ padding: 6px 10px; border-bottom: 1px solid #eee; }}
  tr:last-child td {{ border-bottom: none; }}
  tr.top-row td {{ background: #e6f4f4; }}
  tr:not(.top-row):hover td {{ background: #f5f5f5; }}
  .center {{ text-align: center; }}
  .num    {{ font-variant-numeric: tabular-nums; }}
  .small  {{ font-size: 0.75rem; color: #555; }}
  code    {{ font-family: 'Courier New', monospace; font-size: 0.8rem; }}
  code.seq {{ color: #0d7a7a; font-size: 0.78rem; }}
  .tox-fail {{ color: #c0392b; font-weight: bold; }}
  .tox-warn {{ color: #e67e22; }}
  .ot-hits  {{ color: #c0392b; }}
  .summary  {{ display: flex; gap: 20px; margin-bottom: 20px; flex-wrap: wrap; }}
  .card     {{ background: white; border-radius: 8px; padding: 14px 18px;
               box-shadow: 0 1px 4px rgba(0,0,0,.08); min-width: 160px; }}
  .card-val {{ font-size: 1.5rem; font-weight: bold; color: #0d7a7a; }}
  .card-lab {{ font-size: 0.78rem; color: #888; margin-top: 2px; }}
  .note     {{ font-size: 0.78rem; color: #888; margin-top: 12px; }}
</style>
</head>
<body>
<h1>AlleleSelect v2 <span class="badge">R192Q</span></h1>
<div class="sub">{esc(variant_label)} | ENST00000360228.10 | thexiulab.org</div>

<div class="summary">
  <div class="card">
    <div class="card-val">{len(candidates)}</div>
    <div class="card-lab">Total candidates</div>
  </div>
  <div class="card">
    <div class="card-val">{_fmt(top5[0].get('composite_score','') if top5 else '')}</div>
    <div class="card-lab">Best composite score</div>
  </div>
  <div class="card">
    <div class="card-val">{_fmt(top5[0].get('allele_selectivity_ratio', top5[0].get('allele_selectivity_ratio_kcal_mol','')) if top5 else '')}</div>
    <div class="card-lab">Best ASR (kcal/mol)</div>
  </div>
  <div class="card">
    <div class="card-val">{_fmt(sum(c.get('accessibility_score',0) for c in candidates)/len(candidates) if candidates else 0)}</div>
    <div class="card-lab">Mean accessibility</div>
  </div>
</div>

<table id="tbl">
<thead>
<tr>
  <th onclick="sortTable(0)">Rank</th>
  <th onclick="sortTable(1)">ASO ID</th>
  <th onclick="sortTable(2)">Sequence (5'→3')</th>
  <th onclick="sortTable(3)">Len</th>
  <th onclick="sortTable(4)">ASR (kcal/mol)</th>
  <th onclick="sortTable(5)">Accessibility</th>
  <th onclick="sortTable(6)">Composite</th>
  <th onclick="sortTable(7)">SNP@ASO</th>
  <th onclick="sortTable(8)">Pos score</th>
  <th onclick="sortTable(9)">Region</th>
  <th onclick="sortTable(10)">Off-targets</th>
  <th onclick="sortTable(11)">Splice</th>
  <th onclick="sortTable(12)">Toxicity</th>
  <th>Gapmer pattern</th>
</tr>
</thead>
<tbody>
{''.join(rows_html)}
</tbody>
</table>

<div class="note">
  Highlighted rows (ranks 1–5) = top composite-score candidates. Click column headers to sort.<br>
  <b>ASR</b>: allele selectivity ratio = ΔG<sub>mutant</sub> − ΔG<sub>wildtype</sub> (more negative = better).<br>
  <b>Composite</b> = 0.40×ASR_norm + 0.35×SNP_position_score + 0.25×accessibility.<br>
  <b>SNP@ASO</b>: 1-indexed position of SNP within ASO from 5' end.<br>
  <b>Off-targets</b>: BLASTn hits vs GENCODE v44; "unscreened" = not in top-50 pre-ranking.<br>
  Pipeline v2 | Ostergaard 2013 position scoring | Hagedorn 2022 toxicity screening
</div>

<script>
let sortDir = {{}};
function sortTable(col) {{
  const tbl = document.getElementById('tbl');
  const tbody = tbl.tBodies[0];
  const rows = Array.from(tbody.rows);
  sortDir[col] = !sortDir[col];
  rows.sort((a,b) => {{
    let va = a.cells[col].innerText.trim();
    let vb = b.cells[col].innerText.trim();
    let na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return sortDir[col] ? na-nb : nb-na;
    return sortDir[col] ? va.localeCompare(vb) : vb.localeCompare(va);
  }});
  rows.forEach(r => tbody.appendChild(r));
}}
</script>
</body>
</html>"""

    with open(path, "w", encoding="utf-8") as f:
        f.write(html_out)
    print(f"HTML report saved: {path}")