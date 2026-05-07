"""
alleleselect/scoring/snp_position.py  (v2.0)

SNP Position Scoring + Toxic Sequence Screening
================================================
Based on:
  - Ostergaard et al. 2013 (PMID 23963702): SNP position within DNA gap
    critically determines RNase H allele discrimination in gapmers.
    The mismatch at the gap center maximally disrupts the preferred RNase H
    cleavage site; mismatches at the gap edges or in the wings have far less
    effect on discrimination.
  - Magner et al. 2017 (PMID 28970564): mismatch position and type jointly
    determine SNP-preferential RNase H cleavage.
  - Hagedorn et al. 2022 (PMID 35166597): sequence features predicting acute
    CNS neurotoxicity after ICV dosing. G-rich sequences (especially 3'-end)
    and high PS content are major risk factors.
  - O'Rourke et al. 2026 (PMID 40163297): expanded characterization of acute
    inhibition response; G-rich, C-poor sequences cause the worst outcomes.
  - Burdick et al. 2014: TCC/TGC trinucleotides associated with LNA
    hepatotoxicity.
  - Classic literature: CpG -> immunostimulation (PS backbone); GGGG ->
    G-quadruplex -> non-specific protein binding.
"""

from __future__ import annotations


# ─── SNP POSITION SCORING ─────────────────────────────────────────────────────

def score_snp_position(
    aso_len: int,
    window_start: int,
    snp_cds_pos: int,
    wing_len: int = 5,
) -> dict:
    """
    Calculate where in the ASO the SNP falls and score that position.

    The ASO is the reverse complement of the mRNA window, so:
        snp_in_aso (0-idx) = (aso_len - 1) - (snp_cds_pos - window_start)

    Gap runs from position wing_len to (aso_len - wing_len - 1), 0-indexed.
    Gap center = (gap_start + gap_end) / 2.

    Score:
        SNP in wing        -> 0.0 (RNase H doesn't act in wings)
        SNP at gap edge    -> ~0.08–0.25
        SNP at gap center  -> 1.0 (optimal discrimination)

    Returns dict:
        snp_pos_in_aso  int  (1-indexed, for human readability)
        snp_pos_score   float [0.0, 1.0]
        snp_region      str  'center' | 'mid' | 'edge' | 'wing'
    """
    snp_offset = snp_cds_pos - window_start        # offset from window start
    snp_0idx   = (aso_len - 1) - snp_offset        # position in ASO (0-indexed)

    gap_start  = wing_len
    gap_end    = aso_len - wing_len - 1
    gap_len    = gap_end - gap_start + 1
    gap_center = (gap_start + gap_end) / 2.0

    if snp_0idx < gap_start or snp_0idx > gap_end:
        return {
            "snp_pos_in_aso": snp_0idx + 1,
            "snp_pos_score":  0.0,
            "snp_region":     "wing",
        }

    dist  = abs(snp_0idx - gap_center)
    score = max(0.0, 1.0 - dist / (gap_len / 2.0))

    if score >= 0.80:
        region = "center"
    elif score >= 0.50:
        region = "mid"
    else:
        region = "edge"

    return {
        "snp_pos_in_aso": snp_0idx + 1,
        "snp_pos_score":  round(score, 3),
        "snp_region":     region,
    }


# ─── TOXIC SEQUENCE SCREENING ─────────────────────────────────────────────────

_TOXIC_RULES = [
    # (name, check_fn, severity, reason, citation)
    ("GGGG", lambda s: "GGGG" in s, "serious",
     "G-quadruplex risk: non-specific protein binding / CNS neurotoxicity",
     "Hagedorn 2022 PMID 35166597; O'Rourke 2026 PMID 40163297"),

    ("poly-A≥5", lambda s: "AAAAA" in s, "warning",
     "Poly-A run ≥5: non-specific binding and exonuclease sensitivity",
     "General ASO design guidelines"),

    ("poly-T≥5", lambda s: "TTTTT" in s, "warning",
     "Poly-T run ≥5: potential hairpin formation",
     "General ASO design guidelines"),

    ("poly-C≥5", lambda s: "CCCCC" in s, "warning",
     "Poly-C run ≥5: i-motif formation risk",
     "General ASO design guidelines"),

    ("TCC×2", lambda s: s.count("TCC") >= 2, "warning",
     "TCC trinucleotide ≥2x: hepatotoxicity-associated motif (LNA context)",
     "Burdick 2014"),

    ("TGC×2", lambda s: s.count("TGC") >= 2, "warning",
     "TGC trinucleotide ≥2x: hepatotoxicity-associated motif (LNA context)",
     "Burdick 2014"),

    ("CpG≥3", lambda s: s.count("CG") >= 3, "warning",
     "≥3 CpG motifs: immunostimulatory risk on phosphorothioate backbone",
     "Krieg 1995; Henry 1997"),

    ("CpG×1", lambda s: s.count("CG") >= 1, "note",
     "CpG motif present: note for PS backbone designs",
     "Krieg 1995"),
]


def screen_toxic(aso_seq: str) -> dict:
    """
    Screen ASO DNA sequence (5'->3') for known toxic / problematic motifs.

    Returns dict:
        flags    list of {motif, severity, reason, citation}
        serious  bool  any 'serious' flag present
        warning  bool  any 'serious' or 'warning' flag present
        summary  str   'PASS' | 'WARN:<motif>' | 'FAIL:<motif>'
    """
    seq   = aso_seq.upper().replace("U", "T")
    flags = []
    seen  = set()

    for name, check_fn, severity, reason, citation in _TOXIC_RULES:
        if name in seen:
            continue
        if check_fn(seq):
            flags.append({
                "motif":    name,
                "severity": severity,
                "reason":   reason,
                "citation": citation,
            })
            seen.add(name)

    serious = any(f["severity"] == "serious" for f in flags)
    warning = any(f["severity"] in ("serious", "warning") for f in flags)

    if not flags:
        summary = "PASS"
    elif serious:
        summary = "FAIL:" + ",".join(f["motif"] for f in flags if f["severity"] == "serious")
    else:
        summary = "WARN:" + ",".join(f["motif"] for f in flags if f["severity"] != "note")
        if summary == "WARN:":
            summary = "NOTE:" + ",".join(f["motif"] for f in flags)

    return {
        "flags":   flags,
        "serious": serious,
        "warning": warning,
        "summary": summary,
    }


# ─── COMPOSITE SCORE ──────────────────────────────────────────────────────────

def composite_score(
    asr: float,
    accessibility: float,
    snp_pos_score: float,
    tox_serious: bool = False,
) -> float:
    """
    Combine ASR, accessibility, and SNP position score into one ranking number.

    Candidates with serious toxicity flags are automatically scored 0.0.

    Weights (literature-informed):
        ASR (allele selectivity ratio):  40%
        SNP position score:              35%  (Ostergaard 2013: dominant factor)
        Accessibility:                   25%  (target site availability)

    ASR is normalized: assumes useful range [-2.0, 0.0] kcal/mol.
    More negative ASR = more selective = higher normalized score.

    Returns float in [0.0, 1.0].
    """
    if tox_serious:
        return 0.0

    asr_norm = min(1.0, max(0.0, (-asr) / 2.0))
    acc_norm = min(1.0, max(0.0, float(accessibility)))
    pos_norm = float(snp_pos_score)

    return round(0.40 * asr_norm + 0.35 * pos_norm + 0.25 * acc_norm, 4)


# ─── DEMO / SELF-TEST ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    SNP_POS = 575  # R192Q c.575G>A

    candidates = [
        # id, sequence, window_start, ASR, accessibility
        ("AS_22_565", "CCCTCAGCGTCTGTAGGTCAAA", 565, -0.983, 0.437),
        ("AS_22_559", "GCGTCTGTAGGTCAAACTCCGT", 559, -0.983, 0.428),
        ("AS_21_565", "CCTCAGCGTCTGTAGGTCAAA",  565, -0.983, 0.414),
        ("AS_21_558", "GTCTGTAGGTCAAACTCCGTC",  558, -0.983, 0.413),
        ("AS_18_567", "CTCAGCGTCTGTAGGTCA",      567, -0.983, 0.406),
    ]

    print("AlleleSelect v2 — SNP Position + Toxicity Annotation")
    print("=" * 72)
    print(f"{'ID':<15} {'SNP@ASO':>7} {'Region':>8} {'PosScore':>9} {'Tox':>12} {'Composite':>10}")
    print("-" * 72)

    for cid, seq, wstart, asr, acc in candidates:
        pos = score_snp_position(len(seq), wstart, SNP_POS)
        tox = screen_toxic(seq)
        comp = composite_score(asr, acc, pos["snp_pos_score"], tox["serious"])
        print(f"{cid:<15} {pos['snp_pos_in_aso']:>7} {pos['snp_region']:>8} "
              f"{pos['snp_pos_score']:>9.3f} {tox['summary']:>12} {comp:>10.4f}")

    print()
    print("KEY FINDING for R192Q candidates:")
    print("  AS_21_558: SNP in MOE wing (pos 4) -> score 0.0 -> remove from shortlist")
    print("  AS_21_565: SNP at exact gap center  -> score 1.0 -> promote to rank 1")
