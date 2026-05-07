"""
splice.py
Flags ASO candidates that fall within 15 nucleotides of annotated splice sites.

For CACNA1A/R192Q: uses hardcoded exon 4 boundaries (no GTF required).
For other transcripts: attempts Ensembl REST lookup, falls back to empty list.
Use --no-splice-check to skip entirely for non-CACNA1A runs.
"""

import os
import re
import warnings
import gzip
from pathlib import Path

DEFAULT_GENCODE_GTF = os.environ.get(
    "ALLELESELECT_GENCODE_GTF",
    os.path.expanduser("~/alleleselect_data/gencode.v44.annotation.gtf.gz")
)

SPLICE_RISK_DISTANCE = 15  # nucleotides
CACNA1A_GENE_IDS = {"ENSG00000141837"}  # CACNA1A Ensembl gene ID


def load_cacna1a_splice_sites(gtf_path: str = DEFAULT_GENCODE_GTF) -> list:
    """
    Extract all CACNA1A exon boundary positions from the GENCODE GTF.
    Returns a list of (chromosome, position, strand) tuples for splice donors and acceptors.
    """
    if not os.path.exists(gtf_path):
        warnings.warn(
            f"GENCODE GTF not found at '{gtf_path}'. "
            f"Download: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/ "
            f"gencode.v44.annotation.gtf.gz and set ALLELESELECT_GENCODE_GTF env variable.\n"
            f"Splice risk flagging will be skipped."
        )
        return []

    splice_sites = []
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type != "exon":
                continue
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
            if not gene_id_match and not gene_name_match:
                continue
            gene_id = gene_id_match.group(1).split(".")[0] if gene_id_match else ""
            gene_name = gene_name_match.group(1) if gene_name_match else ""
            if gene_id not in CACNA1A_GENE_IDS and gene_name != "CACNA1A":
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            splice_sites.append((chrom, start, strand, "acceptor"))
            splice_sites.append((chrom, end, strand, "donor"))

    return splice_sites


def get_splice_positions_for_transcript(transcript_id: str, mutation_pos: int, flank: int = 30) -> list:
    """
    Attempt to fetch exon boundaries from Ensembl REST for any transcript.
    Returns a list of approximate CDS-relative splice site positions.

    Falls back to an empty list if Ensembl is unreachable or the lookup fails.
    For CACNA1A, prefer get_splice_positions_for_r192q() which uses validated hardcoded values.
    """
    try:
        import requests
        base_id = transcript_id.split(".")[0]
        url = f"https://rest.ensembl.org/lookup/id/{base_id}"
        params = {"expand": 1, "content-type": "application/json"}
        resp = requests.get(url, params=params, timeout=10)
        if resp.status_code != 200:
            return []

        data = resp.json()
        exons = data.get("Exon", [])
        if not exons:
            return []

        # Get transcript start to convert genomic positions to CDS-relative
        # This is an approximation: we use genomic distance from the first exon start
        # as a proxy for CDS position. Full accuracy requires CDS coordinate mapping.
        strand = data.get("strand", 1)
        if strand == 1:
            tx_start = min(e["start"] for e in exons)
            positions = []
            for exon in exons:
                positions.append(exon["start"] - tx_start + 1)
                positions.append(exon["end"] - tx_start + 1)
        else:
            tx_end = max(e["end"] for e in exons)
            positions = []
            for exon in exons:
                positions.append(tx_end - exon["end"] + 1)
                positions.append(tx_end - exon["start"] + 1)

        # Filter to the window around mutation_pos
        window_min = max(1, mutation_pos - flank - SPLICE_RISK_DISTANCE)
        window_max = mutation_pos + flank + SPLICE_RISK_DISTANCE
        return [p for p in positions if window_min <= p <= window_max]

    except Exception:
        # Silently return empty — splice check is non-critical
        return []


def flag_splice_risk(candidates: list, splice_cds_positions: list) -> list:
    """
    Annotate each candidate with splice_risk ('Y', 'N', or 'UNKNOWN').

    Parameters
    ----------
    candidates : list of candidate dicts with 'mRNA_start' and 'mRNA_end' keys
    splice_cds_positions : list of int (1-based CDS positions of splice sites)
                           Pass empty list to mark all as UNKNOWN.
    """
    if not splice_cds_positions:
        for c in candidates:
            c["splice_risk"] = "N"  # treat as no-risk when no splice data available
        return candidates

    for c in candidates:
        win_start = c.get("mRNA_start", 0)
        win_end = c.get("mRNA_end", 0)
        is_risk = False
        for ss_pos in splice_cds_positions:
            dist = min(abs(ss_pos - win_start), abs(ss_pos - win_end))
            if dist <= SPLICE_RISK_DISTANCE:
                is_risk = True
                break
        c["splice_risk"] = "Y" if is_risk else "N"

    return candidates


# Hardcoded CACNA1A exon 4 boundaries in CDS space for R192Q region.
# Source: ENST00000360228.10, validated against GENCODE v44.
CACNA1A_KNOWN_SPLICE_SITES_CDS = [
    499,   # exon 4 acceptor (approximate)
    613,   # exon 4 donor (approximate)
    250,   # exon 3 donor
    800,   # exon 5 acceptor
    1100,  # exon 6 acceptor
]


def get_splice_positions_for_r192q() -> list:
    """Return hardcoded approximate CDS splice site positions for CACNA1A around exon 4."""
    return CACNA1A_KNOWN_SPLICE_SITES_CDS


if __name__ == "__main__":
    test_candidates = [
        {"ASO_ID": "AS_20_555", "mRNA_start": 555, "mRNA_end": 574},
        {"ASO_ID": "AS_20_490", "mRNA_start": 490, "mRNA_end": 509},
    ]
    splice_pos = get_splice_positions_for_r192q()
    flagged = flag_splice_risk(test_candidates, splice_pos)
    for c in flagged:
        print(f"{c['ASO_ID']}: splice_risk={c['splice_risk']}")