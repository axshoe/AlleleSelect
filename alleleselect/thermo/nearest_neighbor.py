"""
nearest_neighbor.py
SantaLucia 1998 nearest-neighbor thermodynamics for DNA/DNA duplexes.
Implemented from scratch using Table 2 parameters from:
    SantaLucia, J. Jr. (1998). PNAS 95(4):1460-1465.

All enthalpy (dH) values in kcal/mol.
All entropy (dS) values in cal/mol/K.
Gas constant R = 1.987 cal/mol/K.
Temperature for dG calculation: 37°C = 310.15 K.
"""

from math import log

# --------------------------------------------------------------------------
# SantaLucia 1998, Table 2: DNA/DNA nearest-neighbor parameters
# Format: "XY/YX" where XY is the top strand 5'->3', YX is complement
# dH in kcal/mol, dS in cal/mol/K
# --------------------------------------------------------------------------
NN_PARAMS = {
    # Sequence  dH (kcal/mol)  dS (cal/mol/K)
    "AA/TT":   (-7.9,          -22.2),
    "AT/TA":   (-7.2,          -20.4),
    "TA/AT":   (-7.2,          -21.3),
    "CA/GT":   (-8.5,          -22.7),
    "GT/CA":   (-8.4,          -22.4),
    "CT/GA":   (-7.8,          -21.0),
    "GA/CT":   (-8.2,          -22.2),
    "CG/GC":   (-10.6,         -27.2),
    "GC/CG":   (-9.8,          -24.4),
    "GG/CC":   (-8.0,          -19.9),
    "AC/TG":   (-7.8,          -21.0),  # same as CT/GA reversed
    "TC/AG":   (-8.2,          -22.2),  # same as GA/CT reversed
    "AG/TC":   (-7.8,          -21.0),
    "TG/AC":   (-8.5,          -22.7),
    "TT/AA":   (-7.9,          -22.2),
    "CC/GG":   (-8.0,          -19.9),
}

# Initiation parameters (Table 2, SantaLucia 1998)
# With terminal GC: dH = 0.1 kcal/mol, dS = -2.8 cal/mol/K
# With terminal AT: dH = 2.3 kcal/mol, dS =  4.1 cal/mol/K
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3,  4.1)

R_GAS = 1.987   # cal/mol/K
T_BODY = 310.15  # 37°C in Kelvin


def _complement(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(base.upper(), "N")


def _reverse_complement(seq: str) -> str:
    return "".join(_complement(b) for b in reversed(seq))


def _lookup_nn(dinuc: str) -> tuple:
    """
    Look up nearest-neighbor parameters for a 5'->3' dinucleotide pair.
    Returns (dH, dS). Tries both orientations before raising KeyError.
    """
    key1 = dinuc[:2] + "/" + _reverse_complement(dinuc[:2])
    if key1 in NN_PARAMS:
        return NN_PARAMS[key1]
    # Try complement if not found directly
    rc = _reverse_complement(dinuc[:2])
    key2 = rc + "/" + _reverse_complement(rc)
    if key2 in NN_PARAMS:
        return NN_PARAMS[key2]
    raise KeyError(f"No nearest-neighbor parameter for dinucleotide '{dinuc[:2]}'")


def calc_duplex_thermodynamics(aso_seq: str, target_seq: str, c_total_uM: float = 1.0) -> dict:
    """
    Compute duplex thermodynamic parameters using SantaLucia 1998 nearest-neighbor model.

    Parameters
    ----------
    aso_seq : str
        ASO sequence, 5'->3', DNA alphabet.
    target_seq : str
        Target mRNA window (complementary strand), 5'->3', DNA or RNA alphabet.
        Must be same length as aso_seq and fully complementary for perfect-match calculations.
    c_total_uM : float
        Total strand concentration in micromolar. Default = 1.0 uM.

    Returns
    -------
    dict with keys: dH (kcal/mol), dS (cal/mol/K), dG_37 (kcal/mol), Tm_C (°C),
                    aso_seq, target_seq, n_mismatches
    """
    aso_seq = aso_seq.upper().replace("U", "T")
    target_seq = target_seq.upper().replace("U", "T")

    if len(aso_seq) != len(target_seq):
        raise ValueError(
            f"ASO ({len(aso_seq)} nt) and target ({len(target_seq)} nt) must have equal length."
        )

    n = len(aso_seq)
    dH_total = 0.0
    dS_total = 0.0

    # Initiation correction for 5' terminal base of ASO
    first_base = aso_seq[0]
    if first_base in ("G", "C"):
        dH_total += INIT_GC[0]; dS_total += INIT_GC[1]
    else:
        dH_total += INIT_AT[0]; dS_total += INIT_AT[1]

    # Initiation correction for 3' terminal base of ASO
    last_base = aso_seq[-1]
    if last_base in ("G", "C"):
        dH_total += INIT_GC[0]; dS_total += INIT_GC[1]
    else:
        dH_total += INIT_AT[0]; dS_total += INIT_AT[1]

    # Sum nearest-neighbor contributions
    n_mismatches = 0
    for i in range(n - 1):
        b_aso = aso_seq[i: i + 2]
        b_tgt = target_seq[i: i + 2]
        # Check for mismatch at this position
        if _complement(aso_seq[i]) != target_seq[i]:
            n_mismatches += 1
        key = b_aso + "/" + _reverse_complement(b_aso)
        if key in NN_PARAMS:
            dH, dS = NN_PARAMS[key]
        else:
            # Try reverse complement orientation
            rc_key = _reverse_complement(b_aso) + "/" + b_aso
            if rc_key in NN_PARAMS:
                dH, dS = NN_PARAMS[rc_key]
            else:
                # Fall back to average across GC/AT content
                gc = sum(1 for b in b_aso if b in "GC") / 2.0
                dH = -8.0 * gc + -7.9 * (1 - gc)
                dS = -22.0 * gc + -22.2 * (1 - gc)

        dH_total += dH
        dS_total += dS

    # Check last base for mismatch count
    if _complement(aso_seq[-1]) != target_seq[-1]:
        n_mismatches += 1

    # dG at 37°C
    dS_kcal = dS_total / 1000.0  # convert to kcal/mol/K
    dG_37 = dH_total - T_BODY * dS_kcal

    # Melting temperature: Tm = dH / (dS + R * ln(CT/4)) - 273.15
    c_total_mol = c_total_uM * 1e-6
    dS_total_R = dS_total + R_GAS * log(c_total_mol / 4.0)
    if dS_total_R == 0:
        Tm_C = float("nan")
    else:
        Tm_C = (dH_total * 1000.0) / dS_total_R - 273.15

    return {
        "dH": dH_total,
        "dS": dS_total,
        "dG_37": dG_37,
        "Tm_C": Tm_C,
        "aso_seq": aso_seq,
        "target_seq": target_seq,
        "n_mismatches": n_mismatches,
    }


def design_complementary_aso(target_mrna_window: str) -> str:
    """
    Given a target mRNA window (5'->3', RNA or DNA), return the complementary
    ASO sequence (5'->3', DNA).
    """
    return _reverse_complement(target_mrna_window.upper().replace("U", "T"))


if __name__ == "__main__":
    # Validation against SantaLucia 1998 Table 4, Example 1
    # GCATGC / CGTACG: expected dH=-44.1, dS=-118.7, Tm~54.2°C at 50 nM
    aso = "GCATGC"
    tgt = _reverse_complement(aso)
    result = calc_duplex_thermodynamics(aso, tgt, c_total_uM=0.05)
    print(f"dH = {result['dH']:.1f} kcal/mol (expected ~-44.1)")
    print(f"dS = {result['dS']:.1f} cal/mol/K (expected ~-118.7)")
    print(f"Tm = {result['Tm_C']:.1f} °C (expected ~54.2)")
    print(f"dG_37 = {result['dG_37']:.2f} kcal/mol")
