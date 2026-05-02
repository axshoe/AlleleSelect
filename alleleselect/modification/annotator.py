"""
annotator.py
Gapmer modification pattern annotation for ASO candidates.
Applies published Ionis Pharmaceuticals design rules for:
    - Phosphorothioate (PS) backbone positions
    - 2'-MOE (methoxyethyl) or LNA modifications at flanking positions
    - Flagging of toxic PS contexts (PyPy dinucleotides in flanks)

Gapmer architecture: [5' flank (5 nt, modified)] [gap (10 nt, DNA, PS)] [3' flank (5 nt, modified)]
For shorter oligos, proportional design is applied.

References:
    Crooke, S.T., et al. (2021). Antisense technology: An overview and prospectus.
    Nature Reviews Drug Discovery, 20(6), 427-453.
    Geary, R.S., et al. (2015). Pharmacokinetics, biodistribution, and cell uptake
    of antisense oligonucleotides. Advanced Drug Delivery Reviews, 87, 46-51.
"""

# PyPy (pyrimidine-pyrimidine) contexts in flanks associated with increased PS toxicity
# Specifically CC, TT, CT, TC dinucleotides in the 5' flank
TOXIC_PP_DINUCLEOTIDES = {"CC", "TT", "CT", "TC"}

FLANK_SIZE_DEFAULT = 5   # bases at each end with 2'-MOE or LNA modification
GAP_SIZE_DEFAULT = 10    # central DNA gap for RNase H recruitment
MIN_GC_FOR_LNA = 40.0    # % GC below which LNA recommended over 2'-MOE for affinity


def annotate_gapmer(aso_seq: str, use_lna: bool = None) -> dict:
    """
    Annotate a candidate ASO with recommended gapmer modification pattern.

    Parameters
    ----------
    aso_seq : str
        ASO sequence (5'->3', DNA alphabet, 18-22 nt)
    use_lna : bool or None
        If True, recommend LNA at flank positions.
        If None, decide based on GC content (LNA if GC < MIN_GC_FOR_LNA%).

    Returns
    -------
    dict with keys:
        gapmer_pattern (str): annotation string, e.g. "moemoemoe---ddd---moemoe"
        flank_size (int)
        gap_size (int)
        modification_type (str): "MOE" or "LNA"
        ps_toxicity_flag (bool): True if PyPy context detected in 5' flank
        gc_content (float): % GC of aso_seq
        recommended_synthesis_notation (str): compact Ionis-style notation
    """
    aso_seq = aso_seq.upper().replace("U", "T")
    n = len(aso_seq)

    # Compute GC content
    gc_count = sum(1 for b in aso_seq if b in "GC")
    gc_content = 100.0 * gc_count / n if n > 0 else 0.0

    # Determine modification type
    if use_lna is None:
        mod_type = "LNA" if gc_content < MIN_GC_FOR_LNA else "MOE"
    else:
        mod_type = "LNA" if use_lna else "MOE"

    # Compute flank and gap sizes proportionally
    if n >= 20:
        flank = FLANK_SIZE_DEFAULT
        gap = n - 2 * flank
    elif n == 18:
        flank = 4
        gap = 10
    elif n == 19:
        flank = 4
        gap = 11
    else:
        flank = 4
        gap = n - 8

    gap = max(gap, 7)  # minimum gap size for RNase H activity

    # Build pattern string
    mod_label = "L" if mod_type == "LNA" else "M"  # L = LNA, M = MOE
    pattern_parts = (
        [mod_label.lower()] * flank +      # 5' flank
        ["d"] * gap +                       # DNA gap
        [mod_label.lower()] * flank         # 3' flank
    )
    # Pad/trim to sequence length
    while len(pattern_parts) < n:
        pattern_parts.insert(flank, "d")
    pattern_parts = pattern_parts[:n]
    gapmer_pattern = "".join(pattern_parts)

    # Check for PyPy toxicity in 5' flank
    five_prime_flank = aso_seq[:flank]
    ps_toxicity_flag = False
    for i in range(len(five_prime_flank) - 1):
        dinuc = five_prime_flank[i:i+2]
        if dinuc in TOXIC_PP_DINUCLEOTIDES:
            ps_toxicity_flag = True
            break

    # Synthesis notation (compact): e.g. "5MOE-10DNA-5MOE (all-PS backbone)"
    synthesis_notation = (
        f"{flank}{mod_type}-{gap}DNA-{flank}{mod_type} (all-PS backbone)"
    )

    return {
        "gapmer_pattern": gapmer_pattern,
        "flank_size": flank,
        "gap_size": gap,
        "modification_type": mod_type,
        "ps_toxicity_flag": ps_toxicity_flag,
        "gc_content": round(gc_content, 1),
        "recommended_synthesis_notation": synthesis_notation,
    }


def annotate_all_candidates(candidates: list) -> list:
    """
    Apply gapmer annotation to all candidates in-place.
    Adds 'modification' dict and 'recommended_gapmer_pattern' string to each.
    """
    for c in candidates:
        ann = annotate_gapmer(c.get("aso_seq", ""))
        c["modification"] = ann
        c["recommended_gapmer_pattern"] = ann["recommended_synthesis_notation"]
        c["ps_toxicity_flag"] = ann["ps_toxicity_flag"]
    return candidates


if __name__ == "__main__":
    # Test with a representative 20-mer
    test_aso = "GCTTGCTCTCGGTCTTGCCA"
    result = annotate_gapmer(test_aso)
    print(f"Sequence:  {test_aso}")
    print(f"Pattern:   {result['gapmer_pattern']}")
    print(f"Modification: {result['modification_type']}")
    print(f"GC content: {result['gc_content']}%")
    print(f"PS toxicity flag: {result['ps_toxicity_flag']}")
    print(f"Synthesis: {result['recommended_synthesis_notation']}")
