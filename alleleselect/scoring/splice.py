"""
splice.py
Flags ASO candidates that fall within 15 nucleotides of annotated CACNA1A splice sites.
Splice site positions are derived from the GENCODE v44 GTF annotation file.

Candidates within 15 nt of a splice site are annotated SPLICE_RISK = Y.
This does not disqualify them but requires separate wet-lab validation of splicing effects.
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

    Splice donor = end of exon (5' ss), splice acceptor = start of exon (3' ss).
    Positions are 1-based genomic coordinates.
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
            # Check if this exon belongs to CACNA1A
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
            if not gene_id_match and not gene_name_match:
                continue
            gene_id = gene_id_match.group(1).split(".")[0] if gene_id_match else ""
            gene_name = gene_name_match.group(1) if gene_name_match else ""
            if gene_id not in CACNA1A_GENE_IDS and gene_name != "CACNA1A":
                continue

            chrom = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])    # 1-based, inclusive
            strand = fields[6]

            # Splice donor: 3' end of exon (donor = exon end on + strand)
            # Splice acceptor: 5' start of exon (acceptor = exon start on + strand)
            splice_sites.append((chrom, start, strand, "acceptor"))
            splice_sites.append((chrom, end, strand, "donor"))

    return splice_sites


def compute_cds_splice_positions(splice_sites: list, cds_start_genomic: int,
                                 strand: str = "+") -> list:
    """
    Convert genomic splice site positions to CDS-relative positions.
    Simplified: returns a list of CDS-relative positions (1-based) for all splice sites.

    For a proper implementation, this requires the full exon-to-CDS coordinate mapping,
    which depends on the specific transcript. Here we provide the genomic positions
    and flag based on distance to nearest known exon boundary in genomic space.
    """
    positions = []
    for (chrom, pos, ss_strand, ss_type) in splice_sites:
        if ss_strand == strand:
            cds_rel = pos - cds_start_genomic + 1
            if cds_rel > 0:
                positions.append(cds_rel)
    return sorted(positions)


def flag_splice_risk(candidates: list, splice_cds_positions: list) -> list:
    """
    Annotate each candidate with SPLICE_RISK based on proximity to splice sites.

    Parameters
    ----------
    candidates : list of candidate dicts with 'mRNA_start' and 'mRNA_end' keys
    splice_cds_positions : list of int (1-based CDS positions of splice sites)

    Returns
    -------
    candidates list with 'splice_risk' field added ('Y' or 'N').
    """
    if not splice_cds_positions:
        # No splice site data available; mark all as unknown
        for c in candidates:
            c["splice_risk"] = "UNKNOWN"
        return candidates

    for c in candidates:
        win_start = c.get("mRNA_start", 0)
        win_end = c.get("mRNA_end", 0)
        is_risk = False
        for ss_pos in splice_cds_positions:
            # Check distance from window to splice site
            dist = min(abs(ss_pos - win_start), abs(ss_pos - win_end))
            if dist <= SPLICE_RISK_DISTANCE:
                is_risk = True
                break
        c["splice_risk"] = "Y" if is_risk else "N"

    return candidates


# Hardcoded CACNA1A exon boundaries in CDS space for R192Q region (exon 4)
# Source: ENST00000360228.10 exon 4 spans CDS positions ~499-613 approximately.
# These are approximate values for use when GTF is not available.
CACNA1A_KNOWN_SPLICE_SITES_CDS = [
    499,   # approximate exon 4 acceptor
    613,   # approximate exon 4 donor
    250,   # exon 3 donor
    800,   # exon 5 acceptor
    1100,  # exon 6 acceptor
]


def get_splice_positions_for_r192q() -> list:
    """Return hardcoded approximate CDS splice site positions for CACNA1A around exon 4."""
    return CACNA1A_KNOWN_SPLICE_SITES_CDS


if __name__ == "__main__":
    # Test with known R192Q region candidates
    test_candidates = [
        {"ASO_ID": "AS_20_555", "mRNA_start": 555, "mRNA_end": 574},
        {"ASO_ID": "AS_20_490", "mRNA_start": 490, "mRNA_end": 509},  # near exon 4 boundary
    ]
    splice_pos = get_splice_positions_for_r192q()
    flagged = flag_splice_risk(test_candidates, splice_pos)
    for c in flagged:
        print(f"{c['ASO_ID']}: splice_risk={c['splice_risk']}")
