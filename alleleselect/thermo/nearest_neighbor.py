"""
nearest_neighbor.py
Nearest-neighbor thermodynamics for ASO:mRNA duplex stability.

Two parameter sets available:

  (default) Sugimoto 1995 RNA:DNA hybrid parameters
    Sugimoto et al. (1995) Biochemistry 34(35):11211-16. PMID 7545436.
    Measured from 68 RNA:DNA hybrid sequences. Most appropriate for
    MOE PS gapmer (DNA gap) binding to RNA target.
    Recommended by Frank Bennett (Ionis, personal communication 2026).

  (legacy)  SantaLucia 1998 DNA:DNA parameters
    SantaLucia, J. Jr. (1998). PNAS 95(4):1460-1465.
    Original AlleleSelect v1-v6 parameters. Less accurate for RNA targets
    but included for backward compatibility.

Select via calc_duplex_thermodynamics(..., params="sugimoto") [default]
                                       or params="santalucia"

All enthalpy (dH) values in kcal/mol.
All entropy (dS) values in cal/mol/K.
Temperature for dG calculation: 37C = 310.15 K.
"""

from math import log

R_GAS  = 1.987    # cal/mol/K
T_BODY = 310.15   # 37C in Kelvin


# --------------------------------------------------------------------------
# Sugimoto 1995, Table 2: RNA:DNA hybrid nearest-neighbor parameters
# Top strand = RNA (5'->3'), bottom strand = DNA complement (3'->5')
# Format: "rX·dY / rY·dX" where r = RNA base, d = DNA base
# These are written as the RNA dinucleotide over the DNA complement.
# dH in kcal/mol, dS in cal/mol/K
# --------------------------------------------------------------------------
SUGIMOTO_PARAMS = {
    # RNA 5'->3' dinuc : (dH kcal/mol, dS cal/mol/K)
    # From Table 2, Sugimoto 1995 (PMID 7545436)
    "AA": (-11.5, -36.4),
    "AC": (-7.8,  -21.6),
    "AG": (-7.0,  -19.7),
    "AU": (-8.3,  -23.9),
    "CA": (-10.4, -28.4),
    "CC": (-12.8, -31.9),  # note: Sugimoto uses U not T; we convert T->U on RNA strand
    "CG": (-10.4, -26.9),
    "CU": (-9.1,  -23.5),
    "GA": (-8.6,  -22.9),
    "GC": (-8.0,  -17.1),
    "GG": (-9.3,  -23.2),
    "GU": (-5.9,  -12.3),
    "UA": (-7.8,  -23.2),
    "UC": (-5.5,  -13.5),
    "UG": (-9.0,  -26.1),
    "UU": (-6.6,  -18.4),
}

# Sugimoto 1995 initiation parameters
# No terminal correction published; use SantaLucia values as approximation
SUGIMOTO_INIT_GC = (0.1,  -2.8)
SUGIMOTO_INIT_AT = (2.3,   4.1)


# --------------------------------------------------------------------------
# SantaLucia 1998, Table 2: DNA:DNA nearest-neighbor parameters (legacy)
# Format: "XY/YX" where XY is the top strand 5'->3', YX is complement
# --------------------------------------------------------------------------
SANTALUCIA_PARAMS = {
    "AA/TT":  (-7.9,  -22.2),
    "AT/TA":  (-7.2,  -20.4),
    "TA/AT":  (-7.2,  -21.3),
    "CA/GT":  (-8.5,  -22.7),
    "GT/CA":  (-8.4,  -22.4),
    "CT/GA":  (-7.8,  -21.0),
    "GA/CT":  (-8.2,  -22.2),
    "CG/GC":  (-10.6, -27.2),
    "GC/CG":  (-9.8,  -24.4),
    "GG/CC":  (-8.0,  -19.9),
    "AC/TG":  (-7.8,  -21.0),
    "TC/AG":  (-8.2,  -22.2),
    "AG/TC":  (-7.8,  -21.0),
    "TG/AC":  (-8.5,  -22.7),
    "TT/AA":  (-7.9,  -22.2),
    "CC/GG":  (-8.0,  -19.9),
}

SANTALUCIA_INIT_GC = (0.1,  -2.8)
SANTALUCIA_INIT_AT = (2.3,   4.1)


# --------------------------------------------------------------------------
# Peyret 1999 mismatch corrections (applied on top of either parameter set)
# Values approximate the destabilization caused by non-Watson-Crick pairs.
# --------------------------------------------------------------------------
MISMATCH_CORRECTIONS = {
    # (aso_base, target_base): (delta_dH kcal/mol, delta_dS cal/mol/K)
    # Wobble pairs (less destabilizing)
    ("G", "T"): (1.0,   3.4),   # G:T wobble — least destabilizing
    ("T", "G"): (1.0,   3.4),
    ("G", "U"): (1.0,   3.4),
    ("U", "G"): (1.0,   3.4),
    # Purine-purine mismatches
    ("A", "A"): (3.5,  10.5),
    ("G", "A"): (2.8,   8.0),
    ("A", "G"): (2.8,   8.0),
    ("G", "G"): (4.0,  11.5),
    # Pyrimidine-pyrimidine mismatches
    ("C", "T"): (3.5,  10.6),
    ("T", "C"): (3.5,  10.6),
    ("C", "U"): (3.5,  10.6),
    ("U", "C"): (3.5,  10.6),
    ("C", "C"): (4.0,  11.5),
    ("T", "T"): (4.0,  11.5),
    # Purine-pyrimidine non-Watson-Crick
    ("A", "C"): (3.8,  11.2),
    ("C", "A"): (3.8,  11.2),
    ("T", "A"): (0.0,   0.0),   # Watson-Crick — no correction
    ("A", "T"): (0.0,   0.0),
    ("A", "U"): (0.0,   0.0),
    ("U", "A"): (0.0,   0.0),
    ("G", "C"): (0.0,   0.0),
    ("C", "G"): (0.0,   0.0),
}


def _complement_dna(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A"}.get(base.upper(), "N")


def _reverse_complement_dna(seq: str) -> str:
    return "".join(_complement_dna(b) for b in reversed(seq))

# Alias for backward compatibility (tests and other modules import this name)
_reverse_complement = _reverse_complement_dna


def _get_mismatch_correction(aso_base: str, target_base: str) -> tuple:
    """Return (delta_dH, delta_dS) correction for a mismatch pair."""
    key = (aso_base.upper(), target_base.upper().replace("U", "T"))
    return MISMATCH_CORRECTIONS.get(key, (3.0, 9.0))  # fallback: average mismatch


def design_complementary_aso(target_mrna_window: str) -> str:
    """
    Given a target mRNA window (5'->3', RNA or DNA), return the complementary
    ASO sequence (5'->3', DNA). This is the reverse complement of the target.
    """
    return _reverse_complement_dna(target_mrna_window.upper().replace("U", "T"))


def calc_duplex_thermodynamics(
    aso_seq: str,
    target_seq: str,
    c_total_uM: float = 1.0,
    params: str = "sugimoto",
) -> dict:
    """
    Compute duplex thermodynamic parameters for an ASO:RNA target duplex.

    Parameters
    ----------
    aso_seq     : str, ASO sequence 5'->3', DNA alphabet (T not U)
    target_seq  : str, target mRNA window 5'->3', same length as aso_seq
    c_total_uM  : float, total strand concentration in micromolar
    params      : str, "sugimoto" (default, RNA:DNA hybrid, Bennett-recommended)
                       "santalucia" (legacy, DNA:DNA)

    Returns
    -------
    dict: dH (kcal/mol), dS (cal/mol/K), dG_37 (kcal/mol), Tm_C (C),
          aso_seq, target_seq, n_mismatches, params_used
    """
    aso_seq    = aso_seq.upper().replace("U", "T")
    target_seq = target_seq.upper()

    if len(aso_seq) != len(target_seq):
        raise ValueError(
            f"ASO ({len(aso_seq)} nt) and target ({len(target_seq)} nt) must be equal length."
        )

    # Detect if target_seq is the complement strand (same as aso_seq) vs the mRNA strand.
    # score_candidate_window passes target = revcomp(mRNA_window) = aso_seq for perfect match.
    # In that case, the actual mRNA strand is revcomp(target_seq).
    # We normalize: if target_seq == aso_seq (complement strand), flip it back to mRNA strand.
    # More robust: check if target_seq is the revcomp of aso (i.e. the complement strand).
    # If so, the mRNA strand = revcomp(target_seq).
    # We want mRNA_seq (5'->3') for Sugimoto lookup, and for mismatch detection we compare:
    #   Watson-Crick: aso_seq[i] == complement(mRNA_seq[i]) i.e. aso[i]:mRNA[i] form a WC pair
    # Since mRNA and aso are antiparallel, position i of aso pairs with position (n-1-i) of mRNA
    # Equivalently: aso_seq[i] complements mRNA_seq[n-1-i]
    # When target_seq is passed as revcomp(mRNA) = complement strand:
    #   target_seq[i] = complement(mRNA[n-1-i])
    #   So mRNA[n-1-i] = complement(target_seq[i])
    #   Watson-Crick at position i: aso_seq[i] == complement(complement(target_seq[i])) = target_seq[i]
    # So for complement-strand target: match if aso_seq[i] == target_seq[i]
    # For mRNA-strand target: match if aso_seq[i] == complement(mRNA[i]) ... but that's revcomp
    # The safest approach: always convert target to complement strand (antiparallel to ASO)
    # Watson-Crick on complement strand: aso[i] == target_complement[i] (same base, antiparallel)

    # Determine if target_seq is already the complement strand
    # (aso and complement strand have same sequence when aso = revcomp(mRNA))
    # We use a simple heuristic: if more than half the positions have aso[i] == target[i],
    # it's likely the complement strand; otherwise it may be the mRNA strand.
    # But the cleanest approach: always work with mRNA 5'->3' internally.
    # If target looks like the complement strand (aso[i]==target[i] for most positions),
    # convert it to mRNA by taking revcomp.
    same_count = sum(1 for a, t in zip(aso_seq, target_seq) if a == t.replace("U","T"))
    if same_count > len(aso_seq) * 0.5:
        # target_seq is the complement strand; mRNA = revcomp(target_seq)
        mrna_seq = _reverse_complement_dna(target_seq)
        target_is_complement = True
    else:
        # target_seq is the mRNA strand directly
        mrna_seq = target_seq.replace("U", "T")
        target_is_complement = False

    n = len(aso_seq)
    dH_total = 0.0
    dS_total = 0.0
    n_mismatches = 0

    if params == "sugimoto":
        nn_table = SUGIMOTO_PARAMS
        init_gc  = SUGIMOTO_INIT_GC
        init_at  = SUGIMOTO_INIT_AT
        # Sugimoto top strand = RNA 5'->3'. ASO runs 5'->3', mRNA runs antiparallel 3'->5'.
        # mRNA position paired with aso[i] is mrna_seq[n-1-i] (antiparallel).
        # RNA dinucleotide for NN pair (i, i+1) of ASO = mrna positions (n-1-i, n-2-i)
        # read 5'->3' = mrna_seq[n-2-i : n-i] reversed = mrna_seq[n-2-i], mrna_seq[n-1-i]...
        # Actually the mRNA dinucleotide antiparallel to ASO[i:i+2] is mrna_seq[n-2-i:n-i][::-1]
        # = mrna at positions (n-1-i) then (n-2-i) in 3'->5' direction
        # The Sugimoto 5'->3' RNA dinucleotide = complement of aso_dinuc read 5'->3'
        # = _reverse_complement_dna(aso_seq[i:i+2])
        def get_nn(i):
            rna_dinuc = _reverse_complement_dna(aso_seq[i:i+2]).replace("T", "U")
            return nn_table.get(rna_dinuc, (-8.0, -22.0))
    else:
        # SantaLucia DNA:DNA
        nn_table = SANTALUCIA_PARAMS
        init_gc  = SANTALUCIA_INIT_GC
        init_at  = SANTALUCIA_INIT_AT
        def get_nn(i):
            b_aso = aso_seq[i:i+2]
            key1  = b_aso + "/" + _reverse_complement_dna(b_aso)
            if key1 in nn_table:
                return nn_table[key1]
            rc    = _reverse_complement_dna(b_aso)
            key2  = rc + "/" + b_aso
            if key2 in nn_table:
                return nn_table[key2]
            return (-8.0, -22.0)  # fallback

    # Initiation correction for 5' terminal base of ASO
    first_base = aso_seq[0]
    if first_base in ("G", "C"):
        dH_total += init_gc[0]; dS_total += init_gc[1]
    else:
        dH_total += init_at[0]; dS_total += init_at[1]

    # Initiation correction for 3' terminal base of ASO
    last_base = aso_seq[-1]
    if last_base in ("G", "C"):
        dH_total += init_gc[0]; dS_total += init_gc[1]
    else:
        dH_total += init_at[0]; dS_total += init_at[1]

    # Sum nearest-neighbor contributions + mismatch corrections
    for i in range(n - 1):
        dH, dS = get_nn(i)
        dH_total += dH
        dS_total += dS

        # Mismatch detection: ASO[i] pairs with mRNA[n-1-i] (antiparallel)
        # Watson-Crick: complement(aso_seq[i]) == mrna_seq[n-1-i]
        aso_b  = aso_seq[i].upper()
        mrna_b = mrna_seq[n - 1 - i].upper().replace("U", "T")
        if _complement_dna(aso_b) != mrna_b:
            n_mismatches += 1
            mm_dH, mm_dS = _get_mismatch_correction(aso_b, mrna_b)
            dH_total += mm_dH
            dS_total += mm_dS

    # dG at 37C
    dG_37 = dH_total - T_BODY * (dS_total / 1000.0)

    # Melting temperature (Tm)
    c_strand_M = (c_total_uM * 1e-6) / 4.0
    dS_R = dS_total / 1000.0  # convert to kcal/mol/K
    if abs(dS_R) < 1e-10:
        Tm_K = 999.9
    else:
        Tm_K = dH_total / (dS_R + R_GAS / 1000.0 * log(c_strand_M))
    Tm_C = Tm_K - 273.15

    return {
        "dH":           round(dH_total, 3),
        "dS":           round(dS_total, 3),
        "dG_37":        round(dG_37,    4),
        "Tm_C":         round(Tm_C,     2),
        "aso_seq":      aso_seq,
        "target_seq":   target_seq,
        "n_mismatches": n_mismatches,
        "params_used":  params,
    }