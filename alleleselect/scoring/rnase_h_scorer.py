"""
rnase_h_scorer.py
RNase H1 sequence preference scoring for AlleleSelect v8.

Implements the dinucleotide position weight matrix (PWM) from:
    Kielpinski et al. (2017) Nucleic Acids Research 45(22):12932-12944
    doi: 10.1093/nar/gkx1073

Uses the R4b construct dinucleotide model for human RNase H1.
Per Emilio Harris-Mostert (Erasmus MC, personal communication 2026),
this model was adapted for allele-selective gapmer design by appointing
the cleavage site between nucleotides 7 and 8 of a 9-nt window.

The score predicts how efficiently RNase H1 will cleave an RNA:DNA
duplex at a given position based on the surrounding sequence context.
More negative log2 fold-change = more preferred cleavage = higher score.

Usage in AlleleSelect:
    from alleleselect.scoring.rnase_h_scorer import (
        score_rnase_h_cleavage, get_optimal_cleavage_position
    )

    # Score all possible cleavage positions in a gapmer:DNA gap
    cleavage_scores = score_rnase_h_cleavage(target_rna_seq, gap_start, gap_end)

    # Get the position most likely to be cleaved (best for discrimination)
    optimal_pos = get_optimal_cleavage_position(target_rna_seq, gap_start, gap_end)
"""

# --------------------------------------------------------------------------
# R4b dinucleotide model for human RNase H1
# From Kielpinski 2017 Table S1 / Supplementary Excel File008
# Columns 1-9 = positions 1-9 of the 9-nt window
# Cleavage occurs between positions 7 and 8 (Emilio's adaptation)
# Values are log2 fold changes (more negative = more preferred for cleavage)
# --------------------------------------------------------------------------

# R4b dinucleotide model, 9 positions, log2 fold changes
# Keys are RNA dinucleotides (5'->3' on the RNA strand)
# Columns represent positions 1-9 of the 9-nt window
R4B_DINUCLEOTIDE_PWM = {
    #      pos1          pos2          pos3   pos4   pos5          pos6          pos7          pos8          pos9
    "AA": [-0.005433769,  0.0042834876, 0,     0,     0,           -0.009522305,  0.072130987,  0.1267688294, 0.011363718],
    "AC": [-0.003965332, -0.0131270204, 0.03484197, 0, 0,           0.062748291,  0.126789177,  0.1461437685, 0.162159634],
    "AG": [-0.0003775447, 0.0142523519,-0.007480984, 0, 0,          -0.074673607, -0.071779421, -0.0060930033, 0.065834431],
    "AU": [ 0.001430999, -0.0042520103, 0,     0,     0,            0.045505634,  0.039005913,  0.1020504489, 0.090197303],
    "AT": [ 0.001430999, -0.0042520103, 0,     0,     0,            0.045505634,  0.039005913,  0.1020504489, 0.090197303],
    "CA": [ 0.01088227,  -0.0005771547, 0,     0,    -0.011144808,  0.019866681,  0.112734937,  0.1664675358, 0.023867952],
    "CC": [ 0.001501639, -0.0292780898, 0.008212147, 0.018770994, 0.006235185, 0.036813409, 0.103896886, 0.0820432081, 0.139220915],
    "CG": [-0.002819773, -0.0022027671,-0.032077217, 0,  -0.028542785, -0.048494666, -0.036089654, 0.0078608519, 0.042339635],
    "CU": [ 0.01317336,  -0.0021681812, 0,     0,     0.059238522,  0.061439031,  0.126049786,  0.1214791035, 0.081918838],
    "CT": [ 0.01317336,  -0.0021681812, 0,     0,     0.059238522,  0.061439031,  0.126049786,  0.1214791035, 0.081918838],
    "GA": [-0.006397201,  0.0118728021, 0,     0,     0,           -0.024059662, -0.00503588,  -0.0662179231,-0.154121101],
    "GC": [-0.0202702,   -0.0208115596, 0.01434027,-0.009842453, 0, 0.042361668, 0.037128607,  0.0001916755, 0.052317918],
    "GG": [-0.008585581,  0.0046693018, 0.00366422, 0, 0,          -0.106397757, -0.177483795, -0.1673329525,-0.106955224],
    "GU": [-0.004252534, -0.0084865441, 0,     0,     0,            0.038453389, -0.004292496, -0.0978607476,-0.074362022],
    "GT": [-0.004252534, -0.0084865441, 0,     0,     0,            0.038453389, -0.004292496, -0.0978607476,-0.074362022],
    "UA": [ 0.02323227,   0.0120472279, 0,     0,     0,            0.082731471,  0.130425552,  0.1438499459,-0.016044997],
    "TA": [ 0.02323227,   0.0120472279, 0,     0,     0,            0.082731471,  0.130425552,  0.1438499459,-0.016044997],
    "UC": [-4.352105e-06,-0.0073481328, 0.015023827, 0, 0,          0.093607182,  0.141686594,  0.115182071,  0.121599896],
    "TC": [-4.352105e-06,-0.0073481328, 0.015023827, 0, 0,          0.093607182,  0.141686594,  0.115182071,  0.121599896],
    "UG": [ 0.01412205,   0.0129235148,-0.016716362, 0, 0,          -0.002642735, -0.033763896, -0.0249303414,-0.009859235],
    "TG": [ 0.01412205,   0.0129235148,-0.016716362, 0, 0,          -0.002642735, -0.033763896, -0.0249303414,-0.009859235],
    "UU": [ 0.01770772,  -0.0124365496, 0,     0,     0,            0.122115793,  0.120970014,  0.0617152358,-0.005222877],
    "TT": [ 0.01770772,  -0.0124365496, 0,     0,     0,            0.122115793,  0.120970014,  0.0617152358,-0.005222877],
}

# Cleavage position within the 9-nt window: between positions 7 and 8 (0-indexed: after index 6)
CLEAVAGE_POSITION = 7  # 1-indexed: cleavage after position 7


def score_sequence_window(rna_9mer: str) -> float:
    """
    Score a 9-nucleotide RNA sequence window using the R4b dinucleotide PWM.

    Parameters
    ----------
    rna_9mer : str, 9-nucleotide RNA sequence (U or T both accepted)

    Returns
    -------
    float: sum of log2 fold changes for this window (more negative = more
           preferred for RNase H1 cleavage = higher cleavage efficiency)
    """
    seq = rna_9mer.upper().replace("T", "U")
    if len(seq) != 9:
        return 0.0

    total = 0.0
    # Sum dinucleotide contributions at each position
    # Position i (1-indexed) covers the dinucleotide at seq[i-1:i+1]
    for i in range(8):  # positions 1-8 each have a dinucleotide
        dinuc = seq[i:i+2]
        if dinuc in R4B_DINUCLEOTIDE_PWM:
            # Position is i+1 (1-indexed), column index is i
            val = R4B_DINUCLEOTIDE_PWM[dinuc][i]
            total += val * 0.5  # Kielpinski 2017: dinucleotide values multiplied by 0.5

    return total


def score_rnase_h_cleavage(
    mrna_seq: str,
    gap_start_in_mrna: int,
    gap_end_in_mrna: int,
) -> dict:
    """
    Score RNase H1 cleavage efficiency at every position within the DNA gap.

    Slides a 9-nt window across the gap region of the mRNA sequence,
    assigning cleavage to between positions 7 and 8 of each window.

    Parameters
    ----------
    mrna_seq          : str, full mRNA sequence (or sufficient window around gap)
    gap_start_in_mrna : int, 0-indexed start of DNA gap in mRNA coordinates
    gap_end_in_mrna   : int, 0-indexed end of DNA gap in mRNA coordinates (exclusive)

    Returns
    -------
    dict: {cleavage_site_position: score, ...}
          positions are 0-indexed in mRNA coordinates
          score: more negative = more preferred cleavage
    """
    scores = {}
    seq = mrna_seq.upper()

    # Slide 9-nt window; cleavage at position 7 of window (0-indexed: after pos 6)
    # Window must be within the gap region
    window_size = 9
    cleavage_offset = CLEAVAGE_POSITION - 1  # 0-indexed: index 6 of window

    for window_start in range(
        max(0, gap_start_in_mrna - cleavage_offset),
        min(len(seq) - window_size + 1, gap_end_in_mrna - cleavage_offset + 1)
    ):
        window_end = window_start + window_size
        if window_end > len(seq):
            break

        nine_mer = seq[window_start:window_end]
        score = score_sequence_window(nine_mer)

        # Cleavage site in mRNA coordinates
        cleavage_site = window_start + cleavage_offset
        if gap_start_in_mrna <= cleavage_site < gap_end_in_mrna:
            scores[cleavage_site] = score

    return scores


def get_optimal_cleavage_position(
    mrna_seq: str,
    gap_start_in_mrna: int,
    gap_end_in_mrna: int,
) -> tuple:
    """
    Find the position in the DNA gap where RNase H1 most prefers to cleave.

    Returns (optimal_position, score) where optimal_position is 0-indexed
    in mRNA coordinates. More negative score = more efficient cleavage.
    """
    scores = score_rnase_h_cleavage(mrna_seq, gap_start_in_mrna, gap_end_in_mrna)
    if not scores:
        gap_center = (gap_start_in_mrna + gap_end_in_mrna) // 2
        return gap_center, 0.0

    # Most negative score = most preferred cleavage site
    optimal_pos = min(scores, key=lambda p: scores[p])
    return optimal_pos, scores[optimal_pos]


def compute_snp_rnase_h_position_score(
    snp_pos_in_mrna: int,
    mrna_seq: str,
    gap_start_in_mrna: int,
    gap_end_in_mrna: int,
    sigma: float = 2.0,
) -> float:
    """
    Compute SNP position score using RNase H cleavage site prediction.

    Instead of a fixed triangular peak at gap center (Ostergaard 2013),
    this scores the SNP based on its distance from the predicted optimal
    RNase H cleavage site. SNP closer to the cleavage site = higher score.

    This replaces the static triangular scoring in snp_position.py when
    --rnase-h-scoring is enabled.

    Parameters
    ----------
    snp_pos_in_mrna   : int, 0-indexed SNP position in mRNA
    mrna_seq          : str, mRNA sequence
    gap_start_in_mrna : int, 0-indexed gap start
    gap_end_in_mrna   : int, 0-indexed gap end (exclusive)
    sigma             : float, controls width of Gaussian scoring peak

    Returns
    -------
    float: score in [0, 1], 1.0 = SNP at optimal cleavage site
    """
    import math
    optimal_pos, _ = get_optimal_cleavage_position(
        mrna_seq, gap_start_in_mrna, gap_end_in_mrna
    )
    distance = abs(snp_pos_in_mrna - optimal_pos)
    gap_length = gap_end_in_mrna - gap_start_in_mrna
    # Gaussian decay from optimal position
    score = math.exp(-(distance ** 2) / (2 * sigma ** 2))
    return round(score, 4)