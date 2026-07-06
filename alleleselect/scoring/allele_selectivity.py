"""
allele_selectivity.py
Computes allele selectivity ratio (ASR) for each sliding ASO candidate window.

ASR = dG_mutant - dG_wildtype (kcal/mol)
More negative = stronger preferential binding to mutant allele.
Target threshold: ASR < -1.0 kcal/mol (top candidates: ASR < -1.5 kcal/mol)

v7: rna_params argument added to score_candidate_window and generate_candidate_windows.
    Passes through to calc_duplex_thermodynamics to select Sugimoto (default) or
    SantaLucia thermodynamic parameters.
    Mismatch correction is now handled entirely inside calc_duplex_thermodynamics
    (Peyret 1999 corrections applied at the mismatch position). The separate
    apply_mismatch_correction call has been removed to avoid double-counting.
"""

from alleleselect.thermo.nearest_neighbor import (
    calc_duplex_thermodynamics,
    design_complementary_aso,
    _reverse_complement,
)


ASR_THRESHOLD     = -1.0    # kcal/mol minimum for a signal
ASR_TOP_THRESHOLD = -1.5    # kcal/mol for top-ranked candidates
OPTIMAL_POS_MIN   = 8       # 1-indexed from 5' of ASO
OPTIMAL_POS_MAX   = 12
TERMINAL_PENALTY_RANGE = 3  # penalize mutation within this many positions of either end


def score_candidate_window(
    aso_seq: str,
    mutant_window: str,
    wildtype_window: str,
    c_total_uM: float = 1.0,
    rna_params: str = "sugimoto",
) -> dict:
    """
    Compute thermodynamic parameters for one ASO candidate against both alleles.

    Parameters
    ----------
    aso_seq         : str, ASO sequence 5'->3' DNA
    mutant_window   : str, mutant mRNA target window (same length as aso_seq)
    wildtype_window : str, wildtype mRNA target window (same length as aso_seq)
    c_total_uM      : float, strand concentration for Tm calculation
    rna_params      : str, "sugimoto" (default) or "santalucia"
                      Passed to calc_duplex_thermodynamics.

    Returns
    -------
    dict with thermodynamic scores for both alleles and derived ASR.
    """
    aso_seq    = aso_seq.upper().replace("U", "T")
    mut_target = mutant_window.upper().replace("U", "T")
    wt_target  = wildtype_window.upper().replace("U", "T")

    # Pass the complement strand (antiparallel to ASO) as target.
    # calc_duplex_thermodynamics detects this internally and recovers the mRNA strand.
    mut_complement = _reverse_complement(mut_target)
    wt_complement  = _reverse_complement(wt_target)

    # Thermodynamics for both alleles.
    # Mismatch correction (Peyret 1999) is applied inside calc_duplex_thermodynamics.
    thermo_mut = calc_duplex_thermodynamics(
        aso_seq, mut_complement, c_total_uM, params=rna_params
    )
    thermo_wt = calc_duplex_thermodynamics(
        aso_seq, wt_complement, c_total_uM, params=rna_params
    )

    asr = thermo_mut["dG_37"] - thermo_wt["dG_37"]

    # Find mismatch position (where mut and wt windows differ)
    mismatch_positions = [
        i for i in range(len(mut_target))
        if mut_target[i] != wt_target[i]
    ]

    if mismatch_positions:
        mm_pos = mismatch_positions[0]
        # ASO is reverse complement of mutant, so ASO position = (n-1-mm_pos)
        mutation_pos_in_aso = len(aso_seq) - mm_pos  # 1-indexed from 5' of ASO
    else:
        mutation_pos_in_aso = None

    # Score mutation position optimality
    if mutation_pos_in_aso is not None:
        if (mutation_pos_in_aso <= TERMINAL_PENALTY_RANGE or
                mutation_pos_in_aso >= (len(aso_seq) - TERMINAL_PENALTY_RANGE + 1)):
            pos_penalty = 0.5
        else:
            pos_penalty = 0.0
        optimal = OPTIMAL_POS_MIN <= mutation_pos_in_aso <= OPTIMAL_POS_MAX
    else:
        pos_penalty = 0.0
        optimal = False

    return {
        "aso_seq":                   aso_seq,
        "dG_mutant":                 thermo_mut["dG_37"],
        "dG_wildtype":               thermo_wt["dG_37"],
        "Tm_mutant":                 thermo_mut["Tm_C"],
        "Tm_wildtype":               thermo_wt["Tm_C"],
        "dH_mutant":                 thermo_mut["dH"],
        "dS_mutant":                 thermo_mut["dS"],
        "allele_selectivity_ratio":  round(asr, 4),
        "n_mismatches_wt":           thermo_wt["n_mismatches"],
        "mutation_pos_in_aso":       mutation_pos_in_aso,
        "pos_penalty":               pos_penalty,
        "at_optimal_position":       optimal,
        "meets_threshold":           asr < ASR_THRESHOLD,
        "top_candidate":             asr < ASR_TOP_THRESHOLD,
        "params_used":               rna_params,
        # Legacy field names expected by downstream modules
        "delta_G_mutant_kcal_mol":   round(thermo_mut["dG_37"], 4),
        "delta_G_wildtype_kcal_mol": round(thermo_wt["dG_37"], 4),
        "Tm_mutant_C":               round(thermo_mut["Tm_C"], 2),
        "Tm_wildtype_C":             round(thermo_wt["Tm_C"], 2),
    }


def generate_candidate_windows(
    wt_cds: str,
    mut_cds: str,
    mutation_pos: int,
    aso_lengths: list = None,
    flank: int = 30,
    rna_params: str = "sugimoto",
) -> list:
    """
    Slide ASO windows of multiple lengths across the mutation site.

    Parameters
    ----------
    wt_cds       : str, wildtype CDS sequence
    mut_cds      : str, mutant CDS sequence (same length, one substitution)
    mutation_pos : int, 1-based position of the mutation in the CDS
    aso_lengths  : list of int, ASO lengths to try. Default: [18, 19, 20, 21, 22]
    flank        : int, nucleotides to slide on each side of the mutation site
    rna_params   : str, "sugimoto" (default) or "santalucia"

    Returns
    -------
    List of candidate dicts with thermodynamic scores and metadata.
    """
    if aso_lengths is None:
        aso_lengths = [18, 19, 20, 21, 22]

    pos0 = mutation_pos - 1  # convert to 0-indexed
    candidates = []

    for length in aso_lengths:
        for window_start in range(
            max(0, pos0 - flank),
            min(len(wt_cds) - length, pos0 + flank - length + 1) + 1,
        ):
            window_end = window_start + length
            if window_end > len(wt_cds):
                break

            wt_window  = wt_cds[window_start:window_end]
            mut_window = mut_cds[window_start:window_end]

            aso_seq = _reverse_complement(mut_window)

            scores = score_candidate_window(
                aso_seq, mut_window, wt_window,
                rna_params=rna_params,
            )
            scores["length"]    = length
            scores["mRNA_start"] = window_start + 1   # 1-based
            scores["mRNA_end"]   = window_end
            scores["ASO_ID"]     = f"AS_{length}_{window_start + 1}"

            candidates.append(scores)

    return candidates