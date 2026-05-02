"""
allele_selectivity.py
Computes allele selectivity ratio (ASR) for each sliding ASO candidate window.

ASR = dG_mutant - dG_wildtype (kcal/mol)
More negative = stronger preferential binding to mutant allele.
Target threshold: ASR < -1.0 kcal/mol (top candidates: ASR < -1.5 kcal/mol)

Mutation position within ASO:
- Positions 8-12 (of a 20-mer, 1-indexed from 5' of ASO) are optimal.
- Positions 1-3 or last-3 are penalized (terminal mismatches are less discriminating).
"""

from alleleselect.thermo.nearest_neighbor import calc_duplex_thermodynamics, design_complementary_aso
from alleleselect.thermo.mismatch import apply_mismatch_correction, is_mismatch


ASR_THRESHOLD = -1.0     # kcal/mol minimum for a signal
ASR_TOP_THRESHOLD = -1.5 # kcal/mol for top-ranked candidates
OPTIMAL_POS_MIN = 8      # 1-indexed from 5' of ASO
OPTIMAL_POS_MAX = 12
TERMINAL_PENALTY_RANGE = 3  # penalize mutation within this many positions of either end


def score_candidate_window(
    aso_seq: str,
    mutant_window: str,
    wildtype_window: str,
    c_total_uM: float = 1.0,
) -> dict:
    """
    Compute thermodynamic parameters for one ASO candidate against both alleles.

    Parameters
    ----------
    aso_seq : str
        The ASO sequence (5'->3', DNA). Designed to complement the mutant window.
    mutant_window : str
        The mutant mRNA target window (same length as aso_seq, DNA or RNA).
    wildtype_window : str
        The wildtype mRNA target window (same length as aso_seq, DNA or RNA).
    c_total_uM : float
        Strand concentration for Tm calculation.

    Returns
    -------
    dict with thermodynamic scores for both alleles and derived ASR.
    """
    aso_seq = aso_seq.upper().replace("U", "T")
    mut_target = mutant_window.upper().replace("U", "T")
    wt_target = wildtype_window.upper().replace("U", "T")

    # Complement target sequences for duplex calculation
    from alleleselect.thermo.nearest_neighbor import _reverse_complement
    mut_complement = _reverse_complement(mut_target)
    wt_complement = _reverse_complement(wt_target)

    # Perfect-match thermodynamics (ASO vs mutant - designed to be perfect match)
    thermo_mut = calc_duplex_thermodynamics(aso_seq, mut_complement, c_total_uM)

    # Mismatch thermodynamics (ASO vs wildtype - the mutation position is a mismatch)
    thermo_wt_base = calc_duplex_thermodynamics(aso_seq, wt_complement, c_total_uM)

    # Find the mismatch position (where mut and wt sequences differ)
    mismatch_positions = [
        i for i in range(len(mut_target))
        if mut_target[i] != wt_target[i]
    ]

    if mismatch_positions:
        # Apply Peyret correction at the mismatch position
        mm_pos = mismatch_positions[0]
        aso_base_at_mm = aso_seq[len(aso_seq) - 1 - mm_pos]  # ASO is reverse complement
        wt_base = wt_complement[mm_pos]
        thermo_wt = apply_mismatch_correction(thermo_wt_base, aso_base_at_mm, wt_base)
        mutation_pos_in_aso = len(aso_seq) - mm_pos  # 1-indexed from 5' of ASO
    else:
        thermo_wt = thermo_wt_base
        mutation_pos_in_aso = None

    asr = thermo_mut["dG_37"] - thermo_wt["dG_37"]

    # Score mutation position optimality
    pos_penalty = 0.0
    if mutation_pos_in_aso is not None:
        if mutation_pos_in_aso <= TERMINAL_PENALTY_RANGE or \
           mutation_pos_in_aso >= (len(aso_seq) - TERMINAL_PENALTY_RANGE + 1):
            pos_penalty = 0.5  # kcal/mol penalty for terminal position
        optimal = OPTIMAL_POS_MIN <= mutation_pos_in_aso <= OPTIMAL_POS_MAX
    else:
        optimal = False

    return {
        "aso_seq": aso_seq,
        "dG_mutant": thermo_mut["dG_37"],
        "dG_wildtype": thermo_wt["dG_37"],
        "Tm_mutant": thermo_mut["Tm_C"],
        "Tm_wildtype": thermo_wt["Tm_C"],
        "dH_mutant": thermo_mut["dH"],
        "dS_mutant": thermo_mut["dS"],
        "allele_selectivity_ratio": asr,
        "mutation_pos_in_aso": mutation_pos_in_aso,
        "position_optimal": optimal,
        "position_penalty": pos_penalty,
        "adjusted_asr": asr - pos_penalty,  # penalized for terminal positions
        "meets_threshold": asr < ASR_THRESHOLD,
        "top_candidate": asr < ASR_TOP_THRESHOLD and optimal,
    }


def generate_candidate_windows(
    wt_cds: str,
    mut_cds: str,
    mutation_pos: int,
    aso_lengths: list = None,
    flank: int = 30,
) -> list:
    """
    Slide ASO windows of multiple lengths across the mutation site.

    Parameters
    ----------
    wt_cds : str
        Wildtype CDS sequence.
    mut_cds : str
        Mutant CDS sequence (same length, one substitution).
    mutation_pos : int
        1-based position of the mutation in the CDS.
    aso_lengths : list of int
        ASO lengths to try. Default: [18, 19, 20, 21, 22].
    flank : int
        How far to slide on each side of the mutation site (in nucleotides).

    Returns
    -------
    List of candidate dicts with keys: aso_seq, mRNA_start, mRNA_end, length, + scoring fields.
    """
    if aso_lengths is None:
        aso_lengths = [18, 19, 20, 21, 22]

    pos0 = mutation_pos - 1  # 0-indexed
    candidates = []

    for length in aso_lengths:
        # Slide window such that mutation falls at every position within the window
        for window_start in range(
            max(0, pos0 - flank),
            min(len(wt_cds) - length, pos0 + flank - length + 1) + 1,
        ):
            window_end = window_start + length
            if window_end > len(wt_cds):
                break

            wt_window = wt_cds[window_start:window_end]
            mut_window = mut_cds[window_start:window_end]

            # ASO is reverse complement of the mutant target window
            from alleleselect.thermo.nearest_neighbor import _reverse_complement
            aso_seq = _reverse_complement(mut_window)

            scores = score_candidate_window(aso_seq, mut_window, wt_window)
            scores["length"] = length
            scores["mRNA_start"] = window_start + 1  # 1-based
            scores["mRNA_end"] = window_end
            scores["ASO_ID"] = f"AS_{length}_{window_start + 1}"

            candidates.append(scores)

    # Sort by adjusted_asr (most negative first)
    candidates.sort(key=lambda x: x["adjusted_asr"])
    return candidates


if __name__ == "__main__":
    # Quick sanity check with mock sequences
    mock_wt = "AAGACCGAGAGCAAG"
    mock_mut = "AAGACCAAGAGCAAG"  # G>A at position 7
    from alleleselect.thermo.nearest_neighbor import design_complementary_aso, _reverse_complement
    aso = _reverse_complement(mock_mut)
    result = score_candidate_window(aso, mock_mut, mock_wt)
    print(f"dG_mutant: {result['dG_mutant']:.2f} kcal/mol")
    print(f"dG_wildtype: {result['dG_wildtype']:.2f} kcal/mol")
    print(f"ASR: {result['allele_selectivity_ratio']:.2f} kcal/mol")
