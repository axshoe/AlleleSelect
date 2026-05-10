"""
engineered_mismatch.py
Generates ASO variants with one deliberate additional mismatch against the wildtype allele.

Proposed by van Roon-Mom (LUMC) and validated by Ostergaard 2013 Figure 7.
Adding a deliberate mismatch 2-4 positions from the SNP gives the wildtype allele
TWO mismatches (SNP + deliberate) while the mutant allele has only ONE.
This strategy achieved >100-fold selectivity in the HD ASO context.

Activated by: --extra-mismatch flag in cli.py
"""

from __future__ import annotations
from typing import List, Dict, Any

# Nearest-neighbor mismatch penalties (Peyret 1999, DNA/DNA mismatches)
# Approximate delta-delta-G for introducing a mismatch at a position
# These are context-dependent but we use average values for screening
_MISMATCH_PENALTY = {
    # Most destabilizing to least destabilizing (rough ordering)
    "AA": 1.0, "CC": 1.0, "GG": 0.9, "TT": 0.6,
    "AC": 1.0, "CA": 1.0, "GA": 0.8, "AG": 0.8,
    "CT": 0.8, "TC": 0.8, "GC": 0.5, "CG": 0.5,
    "GT": 0.3, "TG": 0.3,  # wobble — least destabilizing
}

# Complement mapping (DNA)
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Most destabilizing base substitutions to introduce (prefer non-wobble)
_DESTABILIZING_SUBS = {
    # For each ASO base, what substitution creates the most destabilizing mismatch
    # with the wildtype mRNA? Avoid creating another wobble (G-T or T-G).
    "A": "C",  # A:C mismatch is very destabilizing
    "C": "A",  # C:A mismatch
    "G": "A",  # G:A mismatch (not T which would be G:T wobble)
    "T": "C",  # T:C mismatch (not A which would create T:A match)
}


def generate_mismatched_candidates(
    top_candidates: List[Dict[str, Any]],
    wt_cds: str,
    snp_cds_pos: int,
    extra_mismatch_type: str = "auto",
    mismatch_offsets: tuple = (-3, -2, 2, 3),  # positions relative to SNP in ASO (0-indexed offset)
    max_to_generate: int = 5,
) -> List[Dict[str, Any]]:
    """
    For each top candidate ASO, generate variants with one deliberate additional
    mismatch against the wildtype allele.

    The extra mismatch is introduced at 2-4 positions from the SNP, ensuring:
    - Mutant mRNA: 1 mismatch (just the SNP) = still binds well
    - Wildtype mRNA: 2 mismatches (SNP + deliberate) = much weaker binding

    Parameters
    ----------
    top_candidates       : list of top candidate dicts (will process first max_to_generate)
    wt_cds               : str, wildtype CDS
    snp_cds_pos          : int, 1-based CDS position of the variant
    extra_mismatch_type  : str, "auto" = most destabilizing; or specific base "A"/"C"/"G"/"T"
    mismatch_offsets     : tuple of int, positions (relative to SNP in ASO, 0-indexed from 5')
                           to try for the extra mismatch. Negative = 5' of SNP, positive = 3'.
    max_to_generate      : int, number of input candidates to process

    Returns
    -------
    list of new candidate dicts with '_em' suffix in ASO_ID, engineered_mismatch=True,
    extra_mismatch_pos, extra_mismatch_base, and updated ASR-equivalent for the
    two-mismatch scenario.
    """
    new_candidates = []
    pool = top_candidates[:max_to_generate]

    for c in pool:
        aso_seq    = c.get("aso_seq", "")
        aso_id     = c.get("ASO_ID", "")
        win_start  = c.get("mRNA_start", 0)
        snp_pos_in = c.get("snp_pos_in_aso", None)  # 1-indexed

        if not aso_seq or snp_pos_in is None:
            continue

        snp_0 = snp_pos_in - 1  # 0-indexed in ASO

        for offset in mismatch_offsets:
            em_pos_0 = snp_0 + offset  # 0-indexed in ASO

            # Must be within the ASO and not at the SNP position
            if em_pos_0 < 0 or em_pos_0 >= len(aso_seq) or em_pos_0 == snp_0:
                continue

            original_base = aso_seq[em_pos_0]

            # Determine mismatch base
            if extra_mismatch_type == "auto":
                em_base = _DESTABILIZING_SUBS.get(original_base, "C")
            else:
                em_base = extra_mismatch_type.upper()

            if em_base == original_base:
                continue  # no change, skip

            # Build new ASO sequence with deliberate mismatch
            new_seq = aso_seq[:em_pos_0] + em_base + aso_seq[em_pos_0 + 1:]

            # The wildtype mRNA base at this position (complement of original ASO base)
            # = original_base complement = what wildtype RNA has at this position
            wt_rna_base = _COMP.get(original_base, "N")
            # The new ASO base vs. wildtype RNA creates a mismatch
            em_mismatch_type = f"{em_base}:{wt_rna_base}"

            # Estimate additional destabilization for wildtype (approximate)
            extra_penalty = _MISMATCH_PENALTY.get(em_mismatch_type[:2], 0.8)

            # New engineered ASR: original ASR minus the extra penalty on wildtype
            # (making wildtype binding even weaker relative to mutant)
            original_asr = c.get("allele_selectivity_ratio", 0.0)
            engineered_asr = original_asr - extra_penalty  # more negative = more selective

            new_c = dict(c)  # copy all fields
            new_c["ASO_ID"]              = f"{aso_id}_em{offset:+d}"
            new_c["aso_seq"]             = new_seq
            new_c["engineered_mismatch"] = True
            new_c["extra_mismatch_pos"]  = em_pos_0 + 1  # 1-indexed in ASO
            new_c["extra_mismatch_base"] = em_base
            new_c["extra_mismatch_type"] = em_mismatch_type
            new_c["extra_penalty_est"]   = round(extra_penalty, 3)
            new_c["allele_selectivity_ratio"] = round(engineered_asr, 3)
            new_c["em_parent_id"]        = aso_id
            new_c["em_offset_from_snp"]  = offset
            # Mark for re-scoring after generation
            new_c["composite_score"]     = 0.0  # will be recomputed in cli.py

            new_candidates.append(new_c)

    return new_candidates


def annotate_em_candidates(
    em_candidates: List[Dict[str, Any]],
    snp_cds_pos: int,
    wing_len: int = 5,
) -> List[Dict[str, Any]]:
    """
    Re-score SNP position and toxicity for engineered-mismatch candidates.
    The extra mismatch changes the sequence, so toxicity must be re-checked.
    SNP position score remains the same (SNP is at same gap position).
    """
    from alleleselect.scoring.snp_position import screen_toxic

    for c in em_candidates:
        # Toxicity re-screen with new sequence
        tox = screen_toxic(c["aso_seq"])
        c["tox_summary"] = tox["summary"]
        c["tox_serious"] = tox["serious"]
        c["tox_warning"] = tox["warning"]
        c["tox_flags"]   = "; ".join(
            f"{f['motif']}: {f['reason']}" for f in tox["flags"]
        ) if tox["flags"] else ""

    return em_candidates
