"""
snp_position.py
SNP position scoring within ASO gap (Ostergaard 2013 framework) + toxic sequence screening
+ chemical modification recommender (v4, Khvorova / Ostergaard feedback).

v4 addition: recommend_gap_modifications()
  Based on Ostergaard 2013 Figures 3 and 7: placing 2S-dT, FRNA, or S-cEt at
  gap positions flanking the SNP suppresses minor RNase H cleavage sites on the
  wildtype duplex and can achieve >100-fold allele discrimination.
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
    Calculate where in the ASO the SNP falls and score it.

    The SNP offset in the ASO is:
        snp_in_aso = (aso_len - 1) - (snp_cds_pos - window_start)
    because the ASO is the reverse complement of the mRNA window.

    Gap runs from index wing_len to (aso_len - wing_len - 1), 0-indexed.
    Center of gap = (gap_start + gap_end) / 2.
    Score = 1.0 - |snp_pos - gap_center| / (gap_len / 2), clamped [0,1].
    SNP in wing -> score = 0.0.

    Returns dict with keys:
      snp_pos_in_aso  (1-indexed, human-readable)
      snp_pos_score   (float 0-1)
      snp_region      (str: "center" | "mid" | "edge" | "wing")
    """
    snp_offset = snp_cds_pos - window_start
    snp_0idx   = (aso_len - 1) - snp_offset

    gap_start  = wing_len
    gap_end    = aso_len - wing_len - 1
    gap_len    = gap_end - gap_start + 1
    gap_center = (gap_start + gap_end) / 2

    if snp_0idx < gap_start or snp_0idx > gap_end:
        return {
            "snp_pos_in_aso": snp_0idx + 1,
            "snp_pos_score":  0.0,
            "snp_region":     "wing",
        }

    dist  = abs(snp_0idx - gap_center)
    score = max(0.0, 1.0 - dist / (gap_len / 2))

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


# ─── CHEMICAL MODIFICATION RECOMMENDER ────────────────────────────────────────

def recommend_gap_modifications(
    aso_seq: str,
    snp_pos_in_aso: int,
    wing_len: int = 5,
) -> dict:
    """
    Recommend chemical modifications at positions flanking the SNP to suppress
    minor RNase H cleavage sites on the wildtype duplex.

    Based on Ostergaard et al. 2013 (PMID 23963702) Figures 3, 5, and 7:
    - 2S-dT at the gap position ONE 5' of the SNP (p-1) is most effective
      when that position is T in the ASO sequence (>48-fold selectivity)
    - FRNA at gap position TWO 5' of the SNP (p-2) achieves >25-fold
    - S-cEt at p-1 achieves good selectivity but may reduce mutant potency

    Parameters
    ----------
    aso_seq        : str, ASO sequence 5'->3'
    snp_pos_in_aso : int, 1-indexed position of SNP in ASO
    wing_len       : int, length of MOE wing (standard: 5)

    Returns
    -------
    dict with:
      recommendations : list of dicts
      synthesis_note  : str (plain-English for synthesis order)
      caveat          : str
    """
    seq   = aso_seq.upper()
    n     = len(seq)
    snp_0 = snp_pos_in_aso - 1  # 0-indexed

    gap_start = wing_len
    gap_end   = n - wing_len - 1

    recommendations = []
    synthesis_lines = []

    pos_p1 = snp_0 - 1  # one 5' of SNP in ASO
    pos_p2 = snp_0 - 2  # two 5' of SNP in ASO

    # Primary recommendation: position p-1
    if gap_start <= pos_p1 <= gap_end:
        base = seq[pos_p1]
        if base == "T":
            recommendations.append({
                "position_1indexed":          pos_p1 + 1,
                "position_in_gap":            pos_p1 - gap_start + 1,
                "aso_base":                   base,
                "modification":               "2S-dT",
                "rationale": (
                    "2S-dT at gap position immediately 5' of SNP. "
                    "Ostergaard 2013 Fig 3A: >48-fold selectivity by blocking "
                    "minor cleavage sites b and c on wildtype duplex."
                ),
                "expected_fold_improvement": ">48-fold (Ostergaard 2013 Fig 3)",
            })
            synthesis_lines.append(
                f"  Pos {pos_p1+1}: T -> 2S-dT  [RECOMMENDED: >48-fold improvement]"
            )
        else:
            recommendations.append({
                "position_1indexed":          pos_p1 + 1,
                "position_in_gap":            pos_p1 - gap_start + 1,
                "aso_base":                   base,
                "modification":               "S-cEt",
                "rationale": (
                    f"S-cEt at gap position immediately 5' of SNP (base={base}, not T). "
                    "Ostergaard 2013 Fig 6A: S-cEt at position 5 gives good selectivity "
                    "without large potency loss. Avoid S-cEt at position 6."
                ),
                "expected_fold_improvement": "good selectivity (Ostergaard 2013 Fig 6)",
            })
            synthesis_lines.append(
                f"  Pos {pos_p1+1}: {base} -> S-cEt  [position 5 of gap, no T available]"
            )

    # Secondary: position p-2, FRNA
    if gap_start <= pos_p2 <= gap_end:
        base2 = seq[pos_p2]
        recommendations.append({
            "position_1indexed":          pos_p2 + 1,
            "position_in_gap":            pos_p2 - gap_start + 1,
            "aso_base":                   base2,
            "modification":               "FRNA",
            "rationale": (
                "FRNA at gap position 2 steps 5' of SNP. "
                "Ostergaard 2013 Fig 5A: >25-fold selectivity standalone. "
                "C3'-endo sugar pucker blocks RNase H at minor cleavage sites. "
                "Can combine with 2S-dT at p-1 for additive effect."
            ),
            "expected_fold_improvement": ">25-fold standalone (Ostergaard 2013 Fig 5)",
        })
        synthesis_lines.append(
            f"  Pos {pos_p2+1}: {base2} -> FRNA  [secondary; >25-fold as standalone]"
        )

    if not recommendations:
        synthesis_note = (
            "SNP is in wing or at gap edge - modification recommendation N/A. "
            "Redesign so SNP falls at gap center."
        )
    else:
        synthesis_note = (
            "Recommended backbone modifications (Ostergaard 2013):\n"
            + "\n".join(synthesis_lines)
            + "\n\nSpecify these during synthesis order. "
            "Start with the primary recommendation; add FRNA secondarily."
        )

    return {
        "recommendations": recommendations,
        "synthesis_note":  synthesis_note,
        "caveat": (
            "Based on Ostergaard 2013 T:G wobble context (HTT/HD). "
            "Applicability to other mismatch types requires experimental confirmation. "
            "Position effect may shift by +/-1 for different mismatch types."
        ),
    }


# ─── TOXIC SEQUENCE SCREENING ─────────────────────────────────────────────────

CPG_MOTIF                  = "CG"
GQUAD_MOTIF                = "GGGG"
HEPATOTOXIC_TRINUCLEOTIDES = ["TCC", "TGC"]
POLYR_RUNS                 = 5
KNOWN_TOXIC_HEXAMERS       = []


def screen_toxic(aso_seq: str) -> dict:
    """
    Screen an ASO DNA sequence for known toxic motifs.

    Returns dict with:
      summary  : str
      serious  : bool
      warning  : bool
      flags    : list of dicts with keys: motif, reason, severity
    """
    seq     = aso_seq.upper().replace("U", "T")
    flags   = []
    serious = False
    warning = False

    cpg_count = seq.count(CPG_MOTIF)
    if cpg_count >= 3:
        flags.append({"motif": CPG_MOTIF,
                      "reason": f">= 3 CpG ({cpg_count}) -> immunostimulatory risk (Krieg 1995)",
                      "severity": "warning"})
        warning = True
    elif cpg_count > 0:
        flags.append({"motif": CPG_MOTIF,
                      "reason": f"{cpg_count} CpG motif(s) - note if PS backbone",
                      "severity": "note"})

    if GQUAD_MOTIF in seq:
        flags.append({"motif": GQUAD_MOTIF,
                      "reason": "GGGG -> G-quadruplex + acute CNS toxicity (Hagedorn 2022)",
                      "severity": "serious"})
        serious = True

    for tri in HEPATOTOXIC_TRINUCLEOTIDES:
        count = seq.count(tri)
        if count >= 2:
            flags.append({"motif": tri,
                          "reason": f"'{tri}' x{count} -> hepatotoxicity-assoc trinucleotide (Burdick 2014)",
                          "severity": "warning"})
            warning = True

    for nuc in "ACGT":
        run = nuc * POLYR_RUNS
        if run in seq:
            flags.append({"motif": run,
                          "reason": f"Poly-{nuc} >={POLYR_RUNS} -> non-specific binding",
                          "severity": "warning"})
            warning = True

    for hex_seq in KNOWN_TOXIC_HEXAMERS:
        if hex_seq in seq:
            flags.append({"motif": hex_seq,
                          "reason": f"Known toxic hexamer '{hex_seq}'",
                          "severity": "serious"})
            serious = True

    if serious:
        summary = "FAIL:" + ",".join(f["motif"] for f in flags if f["severity"] == "serious")
    elif warning:
        summary = "WARN:" + ",".join(f["motif"] for f in flags if f["severity"] == "warning")
    elif flags:
        summary = f"NOTE:CpGx{cpg_count}"
    else:
        summary = "PASS"

    return {"summary": summary, "serious": serious, "warning": warning, "flags": flags}


# ─── COMPOSITE SCORE ──────────────────────────────────────────────────────────

def composite_score(
    asr: float,
    accessibility: float,
    snp_pos_score: float,
    tox_serious: bool = False,
    diff_accessibility: float = 0.0,
) -> float:
    """
    Composite ranking score combining ASR, SNP position, accessibility,
    and (v4) differential accessibility between mutant and wildtype mRNA.

    When diff_accessibility == 0.0 (not computed), uses v2-compatible weights:
      0.40 * ASR_norm + 0.35 * pos_score + 0.25 * accessibility

    When diff_accessibility is provided (v4 mode):
      0.35 * ASR_norm + 0.30 * pos_score + 0.20 * accessibility + 0.15 * diff_norm

    Serious toxicity -> score = 0.0.
    """
    if tox_serious:
        return 0.0

    asr_norm = min(1.0, max(0.0, (-asr) / 2.0))
    acc_norm = min(1.0, max(0.0, accessibility))
    pos_norm = snp_pos_score

    if diff_accessibility == 0.0:
        return round(0.40 * asr_norm + 0.35 * pos_norm + 0.25 * acc_norm, 4)

    # Differential: range roughly -0.5 to +0.5 -> normalize to 0-1
    diff_norm = min(1.0, max(0.0, diff_accessibility + 0.5))
    return round(
        0.35 * asr_norm + 0.30 * pos_norm + 0.20 * acc_norm + 0.15 * diff_norm, 4
    )