"""
offtarget.py
BLASTn-based off-target assessment for ASO candidates against the human transcriptome.
Uses GENCODE v44 human transcriptome FASTA (all protein-coding and lncRNA transcripts).

BLASTn short sequence mode is appropriate for 18-22 nt query sequences.

v4 addition: minimum mismatch distance reporting (Wagner lab feedback).
Reports min_off_target_mismatches and nearest_off_target_gene in addition to hit count.
An ASO with only 1 mismatch to a neuronal off-target is dangerous regardless of hit count.
"""

import csv
import os
import subprocess
import tempfile
import warnings
from pathlib import Path

DEFAULT_BLAST_DB = os.environ.get(
    "ALLELESELECT_BLAST_DB",
    os.path.expanduser("~/alleleselect_data/gencode_v44_transcripts"),
)

OFF_TARGET_IDENTITY   = 80   # % identity threshold
OFF_TARGET_QCOV       = 70   # % query coverage threshold
OFF_TARGET_MIN_LENGTH = 14   # minimum alignment length (nt)

# Risk classification by minimum mismatch count
MISMATCH_RISK_SERIOUS  = 1   # ≤1 mismatch to any off-target = high risk
MISMATCH_RISK_MODERATE = 2   # ≤2 mismatches = moderate risk
MISMATCH_RISK_LOW      = 3   # ≤3 mismatches = low risk (flag but don't disqualify)


def check_blast_available() -> bool:
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def build_blast_db(fasta_path: str, db_path: str = DEFAULT_BLAST_DB) -> None:
    """
    Build a local BLASTn database from a GENCODE v44 transcript FASTA.
    Only needs to be run once. Takes ~5-10 minutes.
    """
    cmd = [
        "makeblastdb",
        "-in", fasta_path,
        "-dbtype", "nucl",
        "-out", db_path,
        "-title", "GENCODE_v44_transcripts",
        "-parse_seqids",
    ]
    print(f"Building BLAST database at {db_path}...")
    subprocess.run(cmd, check=True)
    print("BLAST database built successfully.")


def _compute_mismatches_from_row(row: list) -> int:
    """
    Extract mismatch count from a BLAST tabular output row.
    Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    The 'mismatch' field (index 4) is the raw count of mismatching nucleotides.
    """
    try:
        return int(row[4])
    except (IndexError, ValueError):
        return 999


def run_blast_offtarget(
    candidates: list,
    db_path: str = DEFAULT_BLAST_DB,
    top_n: int = 50,
    gene_name: str = "CACNA1A",
) -> list:
    """
    Run BLASTn short-sequence mode for each top candidate against the human transcriptome.

    v4: Now also computes minimum mismatch distance to any off-target sequence.
    An ASO with 1 mismatch to a highly expressed neuronal gene is high-risk even if
    the overall hit count is low.

    Parameters
    ----------
    candidates : list of dicts (from allele_selectivity.generate_candidate_windows)
    db_path    : str, path to local BLAST database prefix
    top_n      : int, only check the top_n candidates by allele_selectivity_ratio
    gene_name  : str, gene name to exclude from off-target count (self-hits)

    Returns
    -------
    candidates list with these fields added:
      off_target_count          : int, number of off-target hits (excluding self-hits)
      off_target_genes          : list of str, gene names
      min_off_target_mismatches : int, minimum mismatches to any off-target (999 if none)
      nearest_off_target_gene   : str, gene with fewest mismatches ("none" if no hits)
      ot_risk_level             : str, "high" / "moderate" / "low" / "clean" / "unscreened"
    """
    if not check_blast_available():
        warnings.warn(
            "BLAST+ not available. Install from: "
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n"
            "off_target_count will be set to -1 (unknown) for all candidates."
        )
        for c in candidates:
            _set_unscreened(c)
        return candidates

    if not os.path.exists(db_path + ".nhr"):
        warnings.warn(
            f"BLAST database not found at '{db_path}'. "
            f"Run build_blast_db() first. off_target_count set to -1."
        )
        for c in candidates:
            _set_unscreened(c)
        return candidates

    to_blast  = candidates[:top_n]
    remaining = candidates[top_n:]

    work_dir   = tempfile.mkdtemp(prefix="alleleselect_blast_")
    query_fasta = os.path.join(work_dir, "query.fasta")
    blast_out   = os.path.join(work_dir, "blast_results.tsv")

    with open(query_fasta, "w") as f:
        for c in to_blast:
            f.write(f">{c['ASO_ID']}\n{c['aso_seq']}\n")

    cmd = [
        "blastn",
        "-query",         query_fasta,
        "-db",            db_path,
        "-task",          "blastn-short",
        "-perc_identity", str(OFF_TARGET_IDENTITY),
        "-qcov_hsp_perc", str(OFF_TARGET_QCOV),
        "-out",           blast_out,
        "-outfmt",        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-num_threads",   "4",
        "-evalue",        "10",
    ]

    try:
        subprocess.run(cmd, timeout=600, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        warnings.warn(f"BLASTn failed: {e.stderr.decode()[:500]}")
        for c in to_blast + remaining:
            _set_unscreened(c)
        return to_blast + remaining

    # Parse results: collect (gene_name, pident, length, mismatches) per ASO
    hit_map = {}   # ASO_ID -> list of hit dicts
    if os.path.exists(blast_out):
        with open(blast_out, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 12:
                    continue
                qid     = row[0]
                sid     = row[1]
                pident  = float(row[2])
                aln_len = int(row[3])
                n_mm    = _compute_mismatches_from_row(row)

                if pident >= OFF_TARGET_IDENTITY and aln_len >= OFF_TARGET_MIN_LENGTH:
                    gene = _extract_gene_from_sid(sid)
                    if qid not in hit_map:
                        hit_map[qid] = []
                    hit_map[qid].append({
                        "gene":       gene,
                        "sid":        sid,          # full subject ID for self-hit detection
                        "identity":   pident,
                        "length":     aln_len,
                        "mismatches": n_mm,
                    })

    # Annotate candidates
    gene_upper = gene_name.upper()
    for c in to_blast:
        all_hits   = hit_map.get(c["ASO_ID"], [])
        # Exclude self-hits: check gene field AND full subject ID for gene name
        real_hits  = [h for h in all_hits
                      if not _is_self_hit(h["gene"], gene_upper)
                      and not _is_self_hit(h.get("sid", ""), gene_upper)]

        c["off_target_count"]  = len(real_hits)
        c["off_target_genes"]  = [h["gene"] for h in real_hits]

        if real_hits:
            nearest = min(real_hits, key=lambda h: h["mismatches"])
            c["min_off_target_mismatches"] = nearest["mismatches"]
            c["nearest_off_target_gene"]   = nearest["gene"]
            c["ot_risk_level"] = _classify_ot_risk(nearest["mismatches"])
        else:
            c["min_off_target_mismatches"] = 999
            c["nearest_off_target_gene"]   = "none"
            c["ot_risk_level"]             = "clean"

    for c in remaining:
        _set_unscreened(c)

    return to_blast + remaining


def _classify_ot_risk(min_mismatches: int) -> str:
    if min_mismatches <= MISMATCH_RISK_SERIOUS:
        return "high"
    elif min_mismatches <= MISMATCH_RISK_MODERATE:
        return "moderate"
    elif min_mismatches <= MISMATCH_RISK_LOW:
        return "low"
    return "clean"


def _set_unscreened(c: dict) -> None:
    c["off_target_count"]          = -1
    c["off_target_genes"]          = []
    c["min_off_target_mismatches"] = -1
    c["nearest_off_target_gene"]   = "unscreened"
    c["ot_risk_level"]             = "unscreened"


def _is_self_hit(gene_field: str, target_gene_upper: str) -> bool:
    """
    Return True if a BLAST hit is a self-hit to the target gene.
    Checks all pipe-delimited fields, gene name prefixes, and Ensembl gene ID patterns.
    """
    gene_up = gene_field.upper()
    # Direct gene name match or prefix match (e.g. "CACNA1A-201" for "CACNA1A")
    if target_gene_upper in gene_up:
        return True
    # Also check if it looks like an Ensembl transcript of the same gene
    # (handles cases where header parser returned the full subject ID)
    for part in gene_up.split("|"):
        if target_gene_upper in part:
            return True
    return False


def _extract_gene_from_sid(sid: str) -> str:
    """
    Extract gene name from GENCODE transcript ID in BLAST subject ID.
    GENCODE format: ENST00000XXXXX.X|ENSG00000XXXXX.X|HAVANA|GENE_NAME-201|GENE_NAME|...
    Returns the full sid if it can't be parsed, so _is_self_hit can still check all fields.
    """
    parts = sid.split("|")
    if len(parts) >= 5:
        # Prefer field 4 (plain gene name) over field 3 (transcript name like CACNA1A-201)
        return parts[4]
    if len(parts) >= 4:
        return parts[3]
    # Return full sid so downstream filtering can inspect it
    return sid