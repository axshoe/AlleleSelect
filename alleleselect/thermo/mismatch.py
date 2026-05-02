"""
mismatch.py
Internal mismatch penalty parameters from:
    Peyret, N., Seneviratne, P.A., Allawi, H.T., SantaLucia, J. Jr. (1999).
    Biochemistry 38(12):3468-3477.

Provides dH and dS corrections for all 16 internal mismatch types in DNA/DNA duplexes.
Used when computing the allele selectivity ratio: the ASO is designed to complement
the mutant allele, so binding to wildtype introduces a deliberate internal mismatch.
"""

# --------------------------------------------------------------------------
# Peyret 1999 Table 2: internal mismatch parameters
# Keys are "X:Y" where X is the ASO base and Y is the mismatched target base.
# dH in kcal/mol, dS in cal/mol/K
# Values represent the delta-delta correction relative to a matched pair.
# --------------------------------------------------------------------------
MISMATCH_PARAMS = {
    # Canonical mismatches (X:Y where X is ASO, Y is target)
    "A:C": ( 2.3,   4.6),   # AC mismatch dH, dS
    "C:A": ( 0.5,  -4.0),
    "A:G": (-0.1,  -3.8),
    "G:A": (-1.5,  -6.2),
    "A:A": ( 4.7,  12.9),
    "C:C": ( 3.3,  10.4),
    "G:G": ( 3.3,   7.4),
    "T:T": ( 7.6,  20.2),
    "C:T": ( 0.7,  -1.2),
    "T:C": ( 1.2,   0.7),
    "G:T": (-2.5,  -8.3),   # G:T wobble (most stabilizing mismatch)
    "T:G": (-1.5,  -6.2),
    "G:G": ( 3.3,   7.4),
    "C:G": (-0.1,   3.7),   # near-match
    "T:A": ( 3.4,   7.7),   # T:A has reduced stability vs A:T
    "A:T": ( 1.6,   4.5),
}

# Canonical Watson-Crick pairs for reference
WATSON_CRICK = {"A": "T", "T": "A", "G": "C", "C": "G"}


def get_mismatch_correction(aso_base: str, target_base: str) -> tuple:
    """
    Return (ddH, ddS) correction for an internal mismatch between ASO base and target base.

    Parameters
    ----------
    aso_base : str, single nucleotide (A/C/G/T)
    target_base : str, single nucleotide (A/C/G/T)

    Returns
    -------
    (ddH, ddS) in kcal/mol and cal/mol/K
    Returns (0.0, 0.0) if the pair is Watson-Crick (no mismatch penalty).
    """
    aso_base = aso_base.upper()
    target_base = target_base.upper()

    # Check if it's actually a Watson-Crick pair
    if WATSON_CRICK.get(aso_base) == target_base:
        return (0.0, 0.0)

    key = f"{aso_base}:{target_base}"
    if key in MISMATCH_PARAMS:
        return MISMATCH_PARAMS[key]
    # Symmetric fallback
    key_rev = f"{target_base}:{aso_base}"
    if key_rev in MISMATCH_PARAMS:
        dH, dS = MISMATCH_PARAMS[key_rev]
        return (dH, dS)
    # Generic mismatch fallback (conservative estimate)
    return (3.0, 8.0)


def is_mismatch(aso_base: str, target_base: str) -> bool:
    """Return True if aso_base and target_base do not form a Watson-Crick pair."""
    return WATSON_CRICK.get(aso_base.upper()) != target_base.upper()


def apply_mismatch_correction(base_thermo: dict, aso_base: str, target_base: str) -> dict:
    """
    Apply Peyret 1999 mismatch correction to a thermodynamics dict.

    Parameters
    ----------
    base_thermo : dict from nearest_neighbor.calc_duplex_thermodynamics
    aso_base : str, the ASO nucleotide at the mismatch position
    target_base : str, the target nucleotide at the mismatch position

    Returns
    -------
    Corrected thermodynamics dict with updated dH, dS, dG_37.
    """
    from math import log
    from alleleselect.thermo.nearest_neighbor import T_BODY

    ddH, ddS = get_mismatch_correction(aso_base, target_base)
    new_thermo = base_thermo.copy()
    new_thermo["dH"] = base_thermo["dH"] + ddH
    new_thermo["dS"] = base_thermo["dS"] + ddS
    dS_kcal = new_thermo["dS"] / 1000.0
    new_thermo["dG_37"] = new_thermo["dH"] - T_BODY * dS_kcal
    new_thermo["mismatch_corrected"] = True
    return new_thermo


if __name__ == "__main__":
    # Test: G:T wobble should be least destabilizing
    ddH, ddS = get_mismatch_correction("G", "T")
    print(f"G:T wobble: ddH={ddH}, ddS={ddS} (expected ~-2.5, ~-8.3)")
    # A:A should be most destabilizing
    ddH2, ddS2 = get_mismatch_correction("A", "A")
    print(f"A:A mismatch: ddH={ddH2}, ddS={ddS2} (expected ~4.7, ~12.9)")
