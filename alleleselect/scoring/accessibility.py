"""
accessibility.py
RNAfold-based mRNA accessibility scoring for ASO candidate windows.

v4 addition: compute_differential_accessibility()
  Runs RNAfold on BOTH the wildtype and mutant CDS sequences.
  Computes per-candidate: mut_accessibility - wt_accessibility.
  A positive differential means the mutant mRNA is MORE single-stranded
  at the ASO binding window, which correlates with better allele selectivity.
  Motivated by Aguti & Zhou 2024 (PMID 38993932).
"""

import os
import re
import subprocess
import tempfile
import warnings
from typing import Optional


def _neutral_result(n: int) -> dict:
    """Return a neutral (uniform 0.5 accessibility) RNAfold result for n bases."""
    return {
        "n":                 n,
        "mfe_structure":     "." * n,
        "mfe":               0.0,
        "per_base_unpaired": [0.5] * n,
    }


def check_rnafold_available() -> bool:
    try:
        subprocess.run(["RNAfold", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def run_rnafold(sequence: str, work_dir: str = None) -> dict:
    """
    Run RNAfold with partition function to get base-pair probability matrix.

    Parameters
    ----------
    sequence : str
        RNA/DNA sequence (T will be converted to U).
    work_dir : str or None
        Directory for output files. Uses tempdir if None.

    Returns
    -------
    dict with keys:
        mfe_structure      (str)
        mfe                (float)
        per_base_unpaired  (list of float, 1-indexed, length = len(sequence))
    """
    seq_rna = sequence.upper().replace("T", "U")
    n       = len(seq_rna)

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="alleleselect_rnafold_")

    # Write to a named FASTA file. When RNAfold reads from a file with sequence
    # named "query", it creates "query_dp.ps" in the same directory reliably
    # across all platforms. Reading from stdin produces unpredictable filenames
    # on Windows (may be "dot.ps", "stdin_dp.ps", or written to process CWD).
    input_fasta = os.path.join(work_dir, "query.fa")
    with open(input_fasta, "w") as f:
        f.write(f">query\n{seq_rna}\n")

    try:
        result = subprocess.run(
            ["RNAfold", "-p", input_fasta],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=work_dir,
        )
    except subprocess.TimeoutExpired:
        warnings.warn("RNAfold timed out. Returning uniform accessibility = 0.5.")
        return {
            "mfe_structure":     "." * n,
            "mfe":               0.0,
            "per_base_unpaired": [0.5] * n,
        }

    # Parse MFE from stdout: ">query\nSEQUENCE\nSTRUCTURE (MFE)"
    mfe         = 0.0
    mfe_struct  = "." * n
    for line in result.stdout.splitlines():
        m = re.match(r"^([.()\[\]{}]+)\s+\(\s*(-?\d+\.\d+)\s*\)", line)
        if m:
            mfe_struct = m.group(1)
            mfe        = float(m.group(2))
            break

    # Parse base-pair probabilities from _dp.ps file
    per_base_unpaired = _parse_dp_ps(work_dir, n)

    return {
        "mfe_structure":     mfe_struct,
        "mfe":               mfe,
        "per_base_unpaired": per_base_unpaired,
    }


def _parse_dp_ps(work_dir: str, n: int) -> list:
    """
    Parse RNAfold dot-plot PostScript file for base-pair probabilities.
    Returns per-position unpaired probabilities (1 = fully single-stranded).
    """
    # Find the _dp.ps file — now always "query_dp.ps" since we write ">query" FASTA
    # Keep fallbacks for backward compatibility with any cached temp dirs
    dp_candidates = [
        os.path.join(work_dir, "query_dp.ps"),   # primary: our named FASTA approach
        os.path.join(work_dir, "stdin_dp.ps"),
        os.path.join(work_dir, "input_dp.ps"),
        os.path.join(work_dir, "dot.ps"),         # some ViennaRNA versions use this
        os.path.join(work_dir, "sequence_dp.ps"),
    ]
    dp_file = None
    for f in dp_candidates:
        if os.path.exists(f):
            dp_file = f
            break

    if dp_file is None:
        # Search the work_dir for any *_dp.ps file
        try:
            for fname in os.listdir(work_dir):
                if fname.endswith("_dp.ps"):
                    dp_file = os.path.join(work_dir, fname)
                    break
        except OSError:
            pass

    if dp_file is None:
        # List what's actually in work_dir to help diagnose
        try:
            found = os.listdir(work_dir)
            warnings.warn(
                f"RNAfold _dp.ps not found in {work_dir}. "
                f"Files present: {found}. "
                "Returning uniform accessibility = 0.5."
            )
        except OSError:
            warnings.warn("RNAfold _dp.ps not found. Returning uniform accessibility = 0.5.")
        return [0.5] * n

    # ubox entries: i j sqrt(prob) ubox
    paired_prob = {}  # position (1-based) -> sum of paired probabilities
    try:
        with open(dp_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 4 and parts[3] == "ubox":
                    try:
                        i = int(parts[0])
                        j = int(parts[1])
                        sq_prob = float(parts[2])
                        prob = sq_prob ** 2
                        paired_prob[i] = paired_prob.get(i, 0.0) + prob
                        paired_prob[j] = paired_prob.get(j, 0.0) + prob
                    except (ValueError, IndexError):
                        continue
    except OSError:
        return [0.5] * n

    per_base = []
    for i in range(1, n + 1):
        unpaired = max(0.0, min(1.0, 1.0 - paired_prob.get(i, 0.0)))
        per_base.append(unpaired)

    return per_base


def compute_window_accessibility(
    per_base_unpaired: list,
    window_start: int,
    window_end: int,
    cds_start_in_sequence: int = 0,
) -> float:
    """
    Compute mean accessibility for a CDS window in the broader sequence context.

    Parameters
    ----------
    per_base_unpaired    : list of float (from run_rnafold, 1-indexed)
    window_start         : int, 1-based CDS start of window
    window_end           : int, 1-based CDS end of window
    cds_start_in_sequence: int, 0-based offset of CDS start within the RNAfold input sequence

    Returns
    -------
    float: mean accessibility [0, 1]
    """
    seq_start = window_start + cds_start_in_sequence
    seq_end   = window_end   + cds_start_in_sequence

    seq_start = max(1, seq_start)
    seq_end   = min(len(per_base_unpaired), seq_end)

    if seq_start > seq_end:
        return 0.5

    vals = per_base_unpaired[seq_start - 1: seq_end]
    return round(sum(vals) / len(vals), 3) if vals else 0.5


def compute_differential_accessibility(
    wt_cds: str,
    mut_cds: str,
    mutation_pos: int,
    candidates: list,
    flank: int = 200,
    work_dir: Optional[str] = None,
) -> list:
    """
    Compute differential mRNA accessibility: mutant - wildtype at each candidate window.

    Motivated by Aguti & Zhou 2024 (PMID 38993932), which showed that secondary
    structure differences between mutant and wildtype mRNA at the mutation site
    contribute to allele selectivity. Candidates where the mutant is MORE
    single-stranded than wildtype (positive differential) are preferred.

    A positive diff_accessibility means the mutant mRNA is more open to ASO binding,
    while the wildtype may be more structured and harder to bind. This independently
    contributes to allele selectivity beyond thermodynamics alone.

    Parameters
    ----------
    wt_cds      : str, wildtype CDS
    mut_cds     : str, mutant CDS
    mutation_pos: int, 1-based CDS position of the mutation
    candidates  : list of candidate dicts (with mRNA_start, mRNA_end)
    flank       : int, window around mutation for RNAfold (default 200 nt)
    work_dir    : str or None

    Returns
    -------
    candidates list with 'diff_accessibility' field added to each candidate.
    Positive = mutant more accessible. Zero if RNAfold unavailable.
    """
    if not check_rnafold_available():
        warnings.warn(
            "RNAfold not available. diff_accessibility will be 0.0 for all candidates."
        )
        for c in candidates:
            c["wt_accessibility"]   = c.get("accessibility_score", 0.5)
            c["mut_accessibility"]  = c.get("accessibility_score", 0.5)
            c["diff_accessibility"] = 0.0
        return candidates

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="alleleselect_diff_")

    wt_dir  = os.path.join(work_dir, "wt")
    mut_dir = os.path.join(work_dir, "mut")
    os.makedirs(wt_dir,  exist_ok=True)
    os.makedirs(mut_dir, exist_ok=True)

    # Extract sequence window around mutation
    pos0 = mutation_pos - 1
    start = max(0, pos0 - flank)
    end   = min(len(wt_cds), pos0 + flank + 1)

    wt_window  = wt_cds[start:end]
    mut_window = mut_cds[start:end]
    window_cds_start = start + 1  # 1-based CDS start of window

    print("[AlleleSelect] Running differential accessibility (WT vs mutant RNAfold)...")
    wt_result  = run_rnafold(wt_window,  work_dir=wt_dir)
    mut_result = run_rnafold(mut_window, work_dir=mut_dir)

    wt_unpaired  = wt_result["per_base_unpaired"]
    mut_unpaired = mut_result["per_base_unpaired"]

    cds_offset = -(window_cds_start - 1)  # for compute_window_accessibility

    for c in candidates:
        wt_acc  = compute_window_accessibility(
            wt_unpaired, c["mRNA_start"], c["mRNA_end"],
            cds_start_in_sequence=cds_offset
        )
        mut_acc = compute_window_accessibility(
            mut_unpaired, c["mRNA_start"], c["mRNA_end"],
            cds_start_in_sequence=cds_offset
        )
        c["wt_accessibility"]   = wt_acc
        c["mut_accessibility"]  = mut_acc
        c["diff_accessibility"] = round(mut_acc - wt_acc, 3)

    n_pos = sum(1 for c in candidates if c.get("diff_accessibility", 0) > 0)
    print(f"  Differential accessibility complete. "
          f"{n_pos}/{len(candidates)} candidates have positive differential "
          f"(mutant more accessible than wildtype).")

    return candidates