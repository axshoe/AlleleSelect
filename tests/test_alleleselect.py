"""
tests/test_alleleselect.py
Unit tests for AlleleSelect pipeline.
Run with: python -m pytest tests/ -v
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# ============================================================
# 1. HGVS Parser Tests
# ============================================================

from alleleselect.sequence.hgvs_parser import (
    parse_hgvs_coding, validate_against_cds, apply_variant_to_sequence, HGVSParseError
)


def test_hgvs_parse_r192q():
    result = parse_hgvs_coding("c.575G>A")
    assert result["position"] == 575
    assert result["ref"] == "G"
    assert result["alt"] == "A"


def test_hgvs_parse_lowercase():
    result = parse_hgvs_coding("c.575g>a")
    assert result["ref"] == "G"
    assert result["alt"] == "A"


def test_hgvs_parse_invalid():
    with pytest.raises(HGVSParseError):
        parse_hgvs_coding("c.575GA")  # missing arrow


def test_hgvs_parse_same_base():
    with pytest.raises(HGVSParseError):
        parse_hgvs_coding("c.575G>G")


def test_validate_against_cds_correct():
    cds = "A" * 574 + "G" + "A" * 100  # G at position 575
    parsed = {"position": 575, "ref": "G", "alt": "A"}
    validate_against_cds(parsed, cds)  # should not raise


def test_validate_against_cds_mismatch():
    cds = "A" * 574 + "C" + "A" * 100  # C at position 575, not G
    parsed = {"position": 575, "ref": "G", "alt": "A"}
    with pytest.raises(HGVSParseError):
        validate_against_cds(parsed, cds)


def test_apply_variant():
    cds = "AAGACCGAGAGC"
    parsed = {"position": 7, "ref": "G", "alt": "A"}
    mutant = apply_variant_to_sequence(cds, parsed)
    assert mutant == "AAGACCAAGAGC"
    assert len(mutant) == len(cds)


# ============================================================
# 2. Nearest-Neighbor Thermodynamics Tests
# ============================================================

from alleleselect.thermo.nearest_neighbor import calc_duplex_thermodynamics, _reverse_complement


def test_reverse_complement():
    assert _reverse_complement("ATCG") == "CGAT"
    assert _reverse_complement("GCATGC") == "GCATGC"  # palindrome


def test_duplex_simple_homopolymer():
    """AT-rich duplex should have lower Tm than GC-rich."""
    aso_at = "AAAATTTT"
    aso_gc = "GGGGCCCC"
    tgt_at = _reverse_complement(aso_at)
    tgt_gc = _reverse_complement(aso_gc)
    res_at = calc_duplex_thermodynamics(aso_at, tgt_at)
    res_gc = calc_duplex_thermodynamics(aso_gc, tgt_gc)
    assert res_gc["Tm_C"] > res_at["Tm_C"], "GC-rich duplex should have higher Tm"


def test_duplex_length_effect():
    """Longer duplex should have more negative dG."""
    aso_short = "GCATGC"
    aso_long = "GCATGCATGCATGCATGCAT"
    tgt_short = _reverse_complement(aso_short)
    tgt_long = _reverse_complement(aso_long)
    res_short = calc_duplex_thermodynamics(aso_short, tgt_short)
    res_long = calc_duplex_thermodynamics(aso_long, tgt_long)
    assert res_long["dG_37"] < res_short["dG_37"], "Longer duplex should have more negative dG"


def test_duplex_mismatches_weaker():
    """Applying Peyret mismatch correction to A:A should increase dG (weaken binding)."""
    from alleleselect.thermo.mismatch import apply_mismatch_correction
    aso = "GCATGCATGCATGCATGCAT"
    perfect_tgt = _reverse_complement(aso)
    res_perfect = calc_duplex_thermodynamics(aso, perfect_tgt)
    # Apply A:A mismatch correction at an interior position
    res_mm = apply_mismatch_correction(res_perfect, "A", "A")
    assert res_mm["dG_37"] > res_perfect["dG_37"], "A:A mismatch correction should increase dG"


def test_duplex_dg_sign():
    """Spontaneous duplex formation: dG should be negative for stable duplexes."""
    aso = "GCGCGCGCGCGCGCGCGCGC"  # 20-mer, all GC
    tgt = _reverse_complement(aso)
    res = calc_duplex_thermodynamics(aso, tgt)
    assert res["dG_37"] < 0, "Stable GC-rich duplex should have negative dG at 37C"


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        calc_duplex_thermodynamics("AAAA", "TTTTTT")  # different lengths


# ============================================================
# 3. Mismatch Penalty Tests
# ============================================================

from alleleselect.thermo.mismatch import get_mismatch_correction, is_mismatch


def test_watson_crick_no_penalty():
    """Watson-Crick pairs should return zero correction."""
    ddH, ddS = get_mismatch_correction("A", "T")
    assert ddH == 0.0
    assert ddS == 0.0


def test_gt_wobble_least_destabilizing():
    """G:T wobble should have the most negative ddH (least destabilizing mismatch)."""
    ddH_gt, _ = get_mismatch_correction("G", "T")
    ddH_aa, _ = get_mismatch_correction("A", "A")
    assert ddH_gt < ddH_aa, "G:T wobble should be less destabilizing than A:A mismatch"


def test_is_mismatch():
    assert is_mismatch("A", "A") is True
    assert is_mismatch("A", "T") is False
    assert is_mismatch("G", "C") is False
    assert is_mismatch("G", "T") is True


# ============================================================
# 4. Allele Selectivity Tests
# ============================================================

from alleleselect.scoring.allele_selectivity import score_candidate_window


def test_asr_negative_for_mutant_targeting():
    """ASO designed for mutant should prefer mutant (negative ASR)."""
    # R192Q context: mutant has A at position 7, wt has G
    mut_window = "AAGACCAAGAGCAAG"  # A at pos 7
    wt_window  = "AAGACCGAGAGCAAG"  # G at pos 7
    aso = _reverse_complement(mut_window)
    result = score_candidate_window(aso, mut_window, wt_window)
    # ASO complements mutant perfectly, so dG_mutant < dG_wildtype => ASR < 0
    assert result["allele_selectivity_ratio"] < 0, "ASR should be negative for mutant-targeting ASO"


def test_asr_structure_keys():
    """Score result should contain all required keys."""
    mut = "AAGACCAAGAGCAAG"
    wt  = "AAGACCGAGAGCAAG"
    aso = _reverse_complement(mut)
    result = score_candidate_window(aso, mut, wt)
    for key in ["dG_mutant", "dG_wildtype", "allele_selectivity_ratio",
                "Tm_mutant", "Tm_wildtype", "mutation_pos_in_aso",
                "meets_threshold", "top_candidate"]:
        assert key in result, f"Missing key: {key}"


# ============================================================
# 5. Gapmer Annotation Tests
# ============================================================

from alleleselect.modification.annotator import annotate_gapmer


def test_gapmer_pattern_length():
    """Gapmer pattern string should match ASO length."""
    aso = "GCATGCATGCATGCATGCAT"  # 20-mer
    result = annotate_gapmer(aso)
    assert len(result["gapmer_pattern"]) == len(aso)


def test_gapmer_gc_low_uses_lna():
    """Low-GC ASO should recommend LNA."""
    aso_at = "AAAATTTTAAAATTTTAAAA"  # all AT, 0% GC
    result = annotate_gapmer(aso_at)
    assert result["modification_type"] == "LNA"


def test_gapmer_gc_high_uses_moe():
    """High-GC ASO should recommend MOE."""
    aso_gc = "GCGCGCGCGCGCGCGCGCGC"  # all GC, 100% GC
    result = annotate_gapmer(aso_gc)
    assert result["modification_type"] == "MOE"


def test_ps_toxicity_detected():
    """ASO with PyPy context at 5' flank should flag PS toxicity."""
    # CC at 5' end of flank
    aso = "CCAATGCATGCATGCATGCA"
    result = annotate_gapmer(aso)
    assert result["ps_toxicity_flag"] is True


# ============================================================
# 6. Accessibility Module Tests
# ============================================================

from alleleselect.scoring.accessibility import compute_window_accessibility, _neutral_result


def test_neutral_result_structure():
    result = _neutral_result(20)
    assert result["n"] == 20
    assert len(result["per_base_unpaired"]) == 20
    assert all(v == 0.5 for v in result["per_base_unpaired"])


def test_window_accessibility_bounds():
    per_base = [0.8] * 50
    acc = compute_window_accessibility(per_base, 10, 25)
    assert 0.0 <= acc <= 1.0


# ============================================================
# 7. Splice Risk Tests
# ============================================================

from alleleselect.scoring.splice import flag_splice_risk, get_splice_positions_for_r192q


def test_splice_risk_flagged():
    """Window near splice site should be flagged Y."""
    splice_pos = [500]
    candidates = [{"mRNA_start": 495, "mRNA_end": 515, "ASO_ID": "test"}]
    result = flag_splice_risk(candidates, splice_pos)
    assert result[0]["splice_risk"] == "Y"


def test_splice_risk_not_flagged():
    """Window far from splice sites should be N."""
    splice_pos = [500]
    candidates = [{"mRNA_start": 555, "mRNA_end": 575, "ASO_ID": "test"}]
    result = flag_splice_risk(candidates, splice_pos)
    assert result[0]["splice_risk"] == "N"


def test_r192q_splice_positions_populated():
    """Hardcoded R192Q splice positions should be non-empty."""
    pos = get_splice_positions_for_r192q()
    assert len(pos) > 0
    assert all(isinstance(p, int) for p in pos)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
