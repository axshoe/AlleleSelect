"""
offtarget.py
BLASTn-based off-target assessment for ASO candidates against the human transcriptome.
Uses GENCODE v44 human transcriptome FASTA (all protein-coding and lncRNA transcripts).

BLASTn short sequence mode is appropriate for 18-22 nt query sequences.
Off-target threshold: >= 80% identity over >= 14 consecutive nucleotides.

Requirements:
    - BLAST+ installed (blastn, makeblastdb): https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    - GENCODE v44 transcriptome FASTA downloaded and path set in config.
"""

import subprocess
import os
import csv
import warnings
import tempfile
from pathlib import Path

DEFAULT_GENCODE_FASTA = os.environ.get(
    "ALLELESELECT_GENCODE_FASTA",
    os.path.expanduser("~/alleleselect_data/gencode.v44.transcripts.fa")
)
DEFAULT_BLAST_DB = os.environ.get(
    "ALLELESELECT_BLAST_DB",
    os.path.expanduser("~/alleleselect_data/gencode_v44_db")
)

OFF_TARGET_IDENTITY = 80.0  # percent
OFF_TARGET_MIN_LENGTH = 14   # consecutive nt
OFF_TARGET_QCOV = 70.0       # query coverage percent


def check_blast_available() -> bool:
    """Return True if blastn and makeblastdb are on PATH."""
    try:
        r1 = subprocess.run(["blastn", "-version"], capture_output=True, timeout=5)
        r2 = subprocess.run(["makeblastdb", "-version"], capture_output=True, timeout=5)
        return r1.returncode == 0 and r2.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def build_blast_db(fasta_path: str = DEFAULT_GENCODE_FASTA,
                   db_path: str = DEFAULT_BLAST_DB) -> bool:
    """
    Build a local BLASTn database from the GENCODE transcriptome FASTA.
    Only runs if the database files do not already exist.

    Returns True on success, False on failure.
    """
    db_dir = os.path.dirname(db_path)
    os.makedirs(db_dir, exist_ok=True)

    # Check if DB already exists
    if os.path.exists(db_path + ".nhr"):
        return True

    if not os.path.exists(fasta_path):
        warnings.warn(
            f"GENCODE FASTA not found at '{fasta_path}'. "
            f"Download from: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/"
            f"gencode.v44.transcripts.fa.gz and set ALLELESELECT_GENCODE_FASTA env variable."
        )
        return False

    cmd = [
        "makeblastdb",
        "-in", fasta_path,
        "-dbtype", "nucl",
        "-out", db_path,
        "-title", "GENCODE_v44_human_transcriptome",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        warnings.warn(f"makeblastdb failed: {result.stderr[:500]}")
        return False
    return True


def run_blast_offtarget(
    candidates: list,
    db_path: str = DEFAULT_BLAST_DB,
    top_n: int = 50,
) -> list:
    """
    Run BLASTn short-sequence mode for each top candidate against the human transcriptome.

    Parameters
    ----------
    candidates : list of dicts (from allele_selectivity.generate_candidate_windows)
    db_path : str, path to local BLAST database prefix
    top_n : int, only check the top_n candidates by allele_selectivity_ratio

    Returns
    -------
    candidates list with 'off_target_count' and 'off_target_genes' fields added.
    """
    if not check_blast_available():
        warnings.warn(
            "BLAST+ not available. Install from: "
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n"
            "off_target_count will be set to -1 (unknown) for all candidates."
        )
        for c in candidates:
            c["off_target_count"] = -1
            c["off_target_genes"] = []
        return candidates

    if not os.path.exists(db_path + ".nhr"):
        warnings.warn(
            f"BLAST database not found at '{db_path}'. "
            f"Run build_blast_db() first. off_target_count set to -1."
        )
        for c in candidates:
            c["off_target_count"] = -1
            c["off_target_genes"] = []
        return candidates

    # Only BLAST the top candidates
    to_blast = candidates[:top_n]
    remaining = candidates[top_n:]

    work_dir = tempfile.mkdtemp(prefix="alleleselect_blast_")
    query_fasta = os.path.join(work_dir, "query.fasta")
    blast_out = os.path.join(work_dir, "blast_results.tsv")

    # Write multi-sequence query FASTA
    with open(query_fasta, "w") as f:
        for c in to_blast:
            f.write(f">{c['ASO_ID']}\n{c['aso_seq']}\n")

    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_path,
        "-task", "blastn-short",
        "-perc_identity", str(OFF_TARGET_IDENTITY),
        "-qcov_hsp_perc", str(OFF_TARGET_QCOV),
        "-out", blast_out,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-num_threads", "4",
        "-evalue", "10",
    ]

    try:
        subprocess.run(cmd, timeout=600, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        warnings.warn(f"BLASTn failed: {e.stderr.decode()[:500]}")
        for c in to_blast + remaining:
            c["off_target_count"] = -1
            c["off_target_genes"] = []
        return to_blast + remaining

    # Parse results
    hit_map = {}  # ASO_ID -> list of (gene_name, identity, length)
    if os.path.exists(blast_out):
        with open(blast_out, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 12:
                    continue
                qid = row[0]
                sid = row[1]
                pident = float(row[2])
                aln_len = int(row[3])
                if pident >= OFF_TARGET_IDENTITY and aln_len >= OFF_TARGET_MIN_LENGTH:
                    gene_name = _extract_gene_from_sid(sid)
                    if qid not in hit_map:
                        hit_map[qid] = []
                    hit_map[qid].append({"gene": gene_name, "identity": pident, "length": aln_len})

    # Annotate candidates
    for c in to_blast:
        hits = hit_map.get(c["ASO_ID"], [])
        # CACNA1A self-hits are expected and should not be counted as off-targets
        real_hits = [h for h in hits if "CACNA1A" not in h["gene"].upper()]
        c["off_target_count"] = len(real_hits)
        c["off_target_genes"] = [h["gene"] for h in real_hits]

    for c in remaining:
        c["off_target_count"] = -1
        c["off_target_genes"] = []

    return to_blast + remaining


def _extract_gene_from_sid(sid: str) -> str:
    """
    Extract gene name from GENCODE transcript ID in BLAST subject ID.
    GENCODE format: ENST00000XXXXX.X|ENSG00000XXXXX.X|HAVANA|GENE_NAME|TRANSCRIPT_NAME|...
    """
    parts = sid.split("|")
    if len(parts) >= 4:
        return parts[3]
    return sid.split(".")[0]
