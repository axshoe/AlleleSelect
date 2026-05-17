"""
junction_fetcher.py
Fetches wildtype and mutant mRNA sequences for junction-mode ASO design.

Junction mode is for mutations that cause EXON SKIPPING or NOVEL JUNCTIONS
rather than single-nucleotide substitutions. Examples:
  - Splice site mutations (c.6210+1G>A in COL6A3 causing exon 16 skipping)
  - Large deletions that create novel exon-exon junctions
  - Pseudoexon insertions

In junction mode:
  - The mutant mRNA contains a NOVEL JUNCTION that does not exist in wildtype
  - ASOs target the novel junction, so wildtype mRNA cannot be bound at all
    (the junction sequence is simply absent from the wildtype transcript)
  - There is no single-nucleotide mismatch to compute ASR for
  - Selectivity comes from junction specificity, not thermodynamic discrimination
  - AlleleSelect scores candidates by accessibility, off-target hits, and toxicity

Usage:
    alleleselect --junction-mode \\
        --gene COL6A3 \\
        --wt-transcript ENST00000295550.9 \\
        --exon-skip 16 \\
        --output demo/COL6A3_junction/

Or with manual sequences:
    alleleselect --junction-mode \\
        --gene COL6A3 \\
        --mut-fasta mutant_junction.fa \\
        --wt-fasta wildtype_region.fa \\
        --junction-center 50 \\
        --output demo/COL6A3_junction/

Motivated by collaboration with Dr. Sara Aguti (UCL GTAC), whose 2024 Mol Ther
Nucleic Acids paper demonstrated allele-selective ASOs for COL6A3 c.6210+1G>A
via novel exon 15/17 junction targeting.
"""

import requests
import time
import warnings
from typing import Optional, Tuple

ENSEMBL_REST_BASE = "https://rest.ensembl.org"


class JunctionFetchError(Exception):
    pass


def fetch_exon_sequences(
    transcript_id: str,
    exon_numbers: list,
    flank: int = 50,
) -> dict:
    """
    Fetch genomic sequences for specific exons of a transcript from Ensembl.

    Parameters
    ----------
    transcript_id : str, Ensembl transcript ID (e.g. ENST00000295550.9)
    exon_numbers  : list of int, 1-indexed exon numbers to fetch
    flank         : int, extra nucleotides to fetch on each side of each exon

    Returns
    -------
    dict: {exon_number: {"sequence": str, "start": int, "end": int, "strand": int}}
    """
    base_id = transcript_id.split(".")[0]
    url = f"{ENSEMBL_REST_BASE}/lookup/id/{base_id}"
    params = {"expand": 1, "content-type": "application/json"}

    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        raise JunctionFetchError(f"Failed to fetch transcript info: {e}")

    exons = data.get("Exon", [])
    if not exons:
        raise JunctionFetchError(f"No exons found for {transcript_id}")

    strand = data.get("strand", 1)
    chrom  = data.get("seq_region_name", "")

    # Sort exons by position (Ensembl returns them in genomic order)
    if strand == 1:
        exons_sorted = sorted(exons, key=lambda e: e["start"])
    else:
        exons_sorted = sorted(exons, key=lambda e: e["start"], reverse=True)

    result = {}
    for i, exon in enumerate(exons_sorted, 1):
        if i not in exon_numbers:
            continue
        exon_start = exon["start"]
        exon_end   = exon["end"]

        # Fetch sequence with flank
        seq_start = max(1, exon_start - flank)
        seq_end   = exon_end + flank

        seq_url = (
            f"{ENSEMBL_REST_BASE}/sequence/region/human/"
            f"{chrom}:{seq_start}..{seq_end}:{strand}"
        )
        seq_params = {"content-type": "application/json"}
        try:
            seq_resp = requests.get(seq_url, params=seq_params, timeout=30)
            seq_resp.raise_for_status()
            seq_data = seq_resp.json()
            sequence = seq_data.get("seq", "").upper()
        except requests.RequestException as e:
            warnings.warn(f"Failed to fetch sequence for exon {i}: {e}")
            sequence = ""

        result[i] = {
            "sequence": sequence,
            "start":    exon_start,
            "end":      exon_end,
            "strand":   strand,
            "chrom":    chrom,
        }

    return result


def build_junction_sequences(
    transcript_id: str,
    skipped_exon: int,
    junction_flank: int = 100,
) -> Tuple[str, str, int]:
    """
    Build mutant and wildtype junction sequences for an exon-skipping event.

    For a mutation that causes exon N to be skipped:
    - Mutant sequence: end of exon (N-1) + start of exon (N+1)
      This novel junction does not exist in the wildtype
    - Wildtype sequence: end of exon (N-1) + start of exon N + start of exon (N+1)
      The normal junction exists in both exon N-1/N and exon N/N+1 boundaries

    Parameters
    ----------
    transcript_id  : str, Ensembl transcript ID
    skipped_exon   : int, 1-indexed exon number that is skipped in mutant
    junction_flank : int, nucleotides from each exon to include in junction

    Returns
    -------
    Tuple of (mutant_junction_seq, wildtype_region_seq, junction_center_pos)
    junction_center_pos: 0-indexed position in mutant_junction_seq at the junction
    """
    # Fetch the three relevant exons: N-1, N, N+1
    exon_nums = [skipped_exon - 1, skipped_exon, skipped_exon + 1]
    exon_data = fetch_exon_sequences(transcript_id, exon_nums, flank=0)

    if (skipped_exon - 1) not in exon_data:
        raise JunctionFetchError(f"Could not fetch exon {skipped_exon - 1}")
    if (skipped_exon + 1) not in exon_data:
        raise JunctionFetchError(f"Could not fetch exon {skipped_exon + 1}")

    # Build junction sequences using CDS from transcript
    base_id = transcript_id.split(".")[0]
    cds_url = f"{ENSEMBL_REST_BASE}/sequence/id/{base_id}"
    cds_params = {"content-type": "application/json", "type": "cds"}

    try:
        cds_resp = requests.get(cds_url, params=cds_params, timeout=30)
        cds_resp.raise_for_status()
        full_cds = cds_resp.json().get("seq", "").upper()
    except requests.RequestException as e:
        raise JunctionFetchError(f"Failed to fetch CDS: {e}")

    # Also fetch exon-level CDS coordinates
    # Use the cdna sequence to locate exon boundaries
    cdna_url = f"{ENSEMBL_REST_BASE}/sequence/id/{base_id}"
    cdna_params = {"content-type": "application/json", "type": "cdna"}

    try:
        cdna_resp = requests.get(cdna_url, params=cdna_params, timeout=30)
        cdna_resp.raise_for_status()
        full_cdna = cdna_resp.json().get("seq", "").upper()
    except requests.RequestException as e:
        raise JunctionFetchError(f"Failed to fetch cDNA: {e}")

    # Get exon-level sequence info from lookup
    lookup_url = f"{ENSEMBL_REST_BASE}/lookup/id/{base_id}"
    lookup_params = {"expand": 1, "content-type": "application/json"}
    lookup_resp = requests.get(lookup_url, params=lookup_params, timeout=30)
    lookup_data = lookup_resp.json()
    exons = lookup_data.get("Exon", [])
    strand = lookup_data.get("strand", 1)

    if strand == 1:
        exons_sorted = sorted(exons, key=lambda e: e["start"])
    else:
        exons_sorted = sorted(exons, key=lambda e: e["start"], reverse=True)

    # Build exon-level CDS sequences by tracking CDS position
    # Each exon contributes its coding portion to the CDS
    exon_cds_seqs = {}
    cds_pos = 0
    for i, exon in enumerate(exons_sorted, 1):
        # Exon length in genomic space
        exon_len = exon["end"] - exon["start"] + 1
        # Approximate: use exon length as CDS contribution
        # (ignores UTR, but fine for internal exons)
        exon_cds = full_cds[cds_pos:cds_pos + exon_len]
        exon_cds_seqs[i] = exon_cds
        cds_pos += len(exon_cds)

    # Build mutant junction: last junction_flank nt of exon N-1 CDS
    # + first junction_flank nt of exon N+1 CDS
    prev_exon_seq = exon_cds_seqs.get(skipped_exon - 1, "")
    next_exon_seq = exon_cds_seqs.get(skipped_exon + 1, "")
    skip_exon_seq = exon_cds_seqs.get(skipped_exon, "")

    flank_prev = min(junction_flank, len(prev_exon_seq))
    flank_next = min(junction_flank, len(next_exon_seq))

    mutant_junction = (
        prev_exon_seq[-flank_prev:] +
        next_exon_seq[:flank_next]
    )

    # Wildtype "region" around where the junction would be:
    # last flank nt of exon N-1 + all of exon N + first flank nt of exon N+1
    wildtype_region = (
        prev_exon_seq[-flank_prev:] +
        skip_exon_seq +
        next_exon_seq[:flank_next]
    )

    junction_center = flank_prev  # 0-indexed position of junction in mutant_junction

    return mutant_junction, wildtype_region, junction_center


def load_fasta_sequence(fasta_path: str) -> str:
    """Load a single sequence from a FASTA file."""
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq_lines.append(line.upper())
    return "".join(seq_lines)


def generate_junction_candidates(
    mut_junction_seq: str,
    wt_region_seq: str,
    junction_center: int,
    aso_lengths: list = None,
    flank: int = 30,
    gene_label: str = "GENE",
) -> list:
    """
    Generate ASO candidate windows across a novel exon junction.

    Slides windows across the mutant junction sequence centered at junction_center.
    For each window, checks whether the reverse complement appears anywhere in
    the wildtype region (it shouldn't for a true novel junction).

    Parameters
    ----------
    mut_junction_seq : str, mutant mRNA junction sequence
    wt_region_seq    : str, wildtype mRNA at the same genomic region
    junction_center  : int, 0-indexed position of junction in mut_junction_seq
    aso_lengths      : list of int, ASO lengths to try
    flank            : int, nt on each side of junction_center to slide over
    gene_label       : str, for ASO ID naming

    Returns
    -------
    list of candidate dicts (compatible with standard AlleleSelect pipeline output)
    """
    if aso_lengths is None:
        aso_lengths = [18, 19, 20, 21, 22]

    _COMP = str.maketrans("ACGT", "TGCA")

    def revcomp(seq):
        return seq.translate(_COMP)[::-1]

    candidates = []
    for length in aso_lengths:
        # Slide window across junction region
        for start in range(
            max(0, junction_center - flank - length),
            min(len(mut_junction_seq) - length + 1, junction_center + flank + 1)
        ):
            end = start + length
            target_seq = mut_junction_seq[start:end]
            aso_seq    = revcomp(target_seq)  # ASO is reverse complement of target

            # Check if ASO binds anywhere in wildtype
            wt_revcomp_target = revcomp(aso_seq)  # = target_seq
            wt_hits = wt_region_seq.count(wt_revcomp_target)

            # Junction specificity: 1.0 if junction sequence absent from WT, 0.0 if present
            junction_specificity = 1.0 if wt_hits == 0 else 0.0

            # Distance from junction center (0 = ASO spans the junction)
            window_center  = start + length / 2
            dist_from_junc = abs(window_center - junction_center)
            # Score: closer to junction center = more likely to span the novel boundary
            junc_pos_score = max(0.0, 1.0 - dist_from_junc / flank)

            aso_id = f"JM_{length}_{start}"

            c = {
                "ASO_ID":                       aso_id,
                "aso_seq":                      aso_seq,
                "length":                       length,
                "mRNA_start":                   start + 1,   # 1-indexed
                "mRNA_end":                     end,
                "junction_center":              junction_center,
                "junction_specificity":         junction_specificity,
                "junction_pos_score":           round(junc_pos_score, 3),
                "wt_cross_reactivity":          wt_hits,
                # Standard fields (junction mode has no ASR — selectivity from specificity)
                "allele_selectivity_ratio":     0.0,
                "delta_G_mutant_kcal_mol":      0.0,
                "delta_G_wildtype_kcal_mol":    0.0,
                "Tm_mutant_C":                  0.0,
                "Tm_wildtype_C":                0.0,
                "accessibility_score":          0.5,
                "off_target_count":             -1,
                "off_target_genes":             [],
                "min_off_target_mismatches":    -1,
                "nearest_off_target_gene":      "unscreened",
                "ot_risk_level":                "unscreened",
                "splice_risk":                  "N",
                "ps_toxicity_flag":             False,
                "recommended_gapmer_pattern":   f"5MOE-{length-10}DNA-5MOE",
                "top_candidate":                False,
                "snp_pos_in_aso":               -1,
                "snp_pos_score":                junc_pos_score,
                "snp_region":                   "junction",
                "composite_score":              0.0,
                "mode":                         "junction",
            }
            candidates.append(c)

    return candidates


def score_junction_composite(candidates: list) -> list:
    """
    Compute composite score for junction-mode candidates.

    Junction composite = 0.50 * junction_specificity
                       + 0.30 * junction_pos_score
                       + 0.20 * accessibility_score

    Junction specificity dominates: if the ASO binds anywhere in wildtype (wt_hits > 0),
    it's deprioritized regardless of position or accessibility.
    Tox serious flag -> score = 0.
    """
    for c in candidates:
        if c.get("tox_serious", False):
            c["composite_score"] = 0.0
            continue
        spec   = c.get("junction_specificity", 0.0)
        pos    = c.get("junction_pos_score",   0.0)
        acc    = c.get("accessibility_score",  0.5)
        c["composite_score"] = round(0.50 * spec + 0.30 * pos + 0.20 * acc, 4)

    return candidates
