"""
accessibility.py
mRNA secondary structure accessibility scoring using RNAfold (Vienna RNA package).
Calls RNAfold via subprocess, parses the base-pair probability matrix (.dp file),
and computes per-window accessibility as 1 - mean_base_pair_probability.

RNAfold reference:
    Lorenz, R., et al. (2011). ViennaRNA Package 2.0. Algorithms Mol Biol 6(1):26.
    Free download: https://www.tbi.univie.ac.at/RNA/
"""

import subprocess
import os
import re
import tempfile
import warnings


def check_rnafold_available() -> bool:
    """Return True if RNAfold is on PATH, False otherwise."""
    try:
        result = subprocess.run(
            ["RNAfold", "--version"],
            capture_output=True, text=True, timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def run_rnafold(sequence: str, work_dir: str = None) -> dict:
    """
    Run RNAfold with partition function to get base-pair probability matrix.

    Parameters
    ----------
    sequence : str
        RNA sequence (A/C/G/U or DNA with T converted to U).
    work_dir : str or None
        Directory for output files. Uses tempdir if None.

    Returns
    -------
    dict with keys:
        mfe_structure (str): MFE dot-bracket structure
        mfe (float): MFE in kcal/mol
        bp_matrix (dict): {(i,j): probability} base-pair probability matrix
        per_base_unpaired (list): per-position unpaired probability (1 - sum of paired probs)
    """
    seq_rna = sequence.upper().replace("T", "U")
    n = len(seq_rna)

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="alleleselect_rnafold_")

    # Write input to temp file
    input_path = os.path.join(work_dir, "input.fasta")
    with open(input_path, "w") as f:
        f.write(f">target\n{seq_rna}\n")

    # Run RNAfold with partition function
    cmd = [
        "RNAfold",
        "-p",           # partition function + base-pair probability matrix
        "--noPS",       # skip PostScript output
        input_path,
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=work_dir,
        )
    except FileNotFoundError:
        warnings.warn(
            "RNAfold not found. Install ViennaRNA: https://www.tbi.univie.ac.at/RNA/\n"
            "Accessibility scores will be set to 0.5 (neutral) for all windows."
        )
        return _neutral_result(n)
    except subprocess.TimeoutExpired:
        warnings.warn("RNAfold timed out. Using neutral accessibility scores.")
        return _neutral_result(n)

    if result.returncode != 0:
        warnings.warn(f"RNAfold error: {result.stderr[:500]}")
        return _neutral_result(n)

    # Parse MFE structure and energy from stdout
    lines = result.stdout.strip().split("\n")
    mfe_structure = ""
    mfe = 0.0
    for line in lines:
        m = re.match(r"^([.()\[\]{}|,]+)\s+\(\s*(-?\d+\.\d+)\)", line)
        if m:
            mfe_structure = m.group(1)
            mfe = float(m.group(2))
            break

    # Parse base-pair probability matrix from _dp.ps file
    dp_file = os.path.join(work_dir, "target_dp.ps")
    bp_matrix = {}
    if os.path.exists(dp_file):
        bp_matrix = _parse_dp_file(dp_file)
    else:
        # Try current directory
        dp_file_alt = "target_dp.ps"
        if os.path.exists(dp_file_alt):
            bp_matrix = _parse_dp_file(dp_file_alt)

    # Compute per-base unpaired probability
    per_base_paired = [0.0] * (n + 1)  # 1-indexed
    for (i, j), prob in bp_matrix.items():
        per_base_paired[i] += prob ** 2  # RNAfold dp.ps stores sqrt of probability
        per_base_paired[j] += prob ** 2

    per_base_unpaired = [
        max(0.0, min(1.0, 1.0 - per_base_paired[i]))
        for i in range(1, n + 1)
    ]

    return {
        "mfe_structure": mfe_structure,
        "mfe": mfe,
        "bp_matrix": bp_matrix,
        "per_base_unpaired": per_base_unpaired,
        "n": n,
    }


def _parse_dp_file(dp_path: str) -> dict:
    """
    Parse RNAfold _dp.ps PostScript file to extract base-pair probability matrix.
    Returns dict: {(i, j): sqrt_probability} where i < j (1-indexed positions).
    """
    bp_matrix = {}
    pattern = re.compile(r"(\d+)\s+(\d+)\s+([\d.e+-]+)\s+ubox")
    try:
        with open(dp_path, "r") as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    i, j, sq_prob = int(m.group(1)), int(m.group(2)), float(m.group(3))
                    bp_matrix[(i, j)] = sq_prob  # stored as sqrt(prob) in dp.ps
    except Exception:
        pass
    return bp_matrix


def _neutral_result(n: int) -> dict:
    """Fallback result when RNAfold is unavailable: all positions 0.5 accessible."""
    return {
        "mfe_structure": "." * n,
        "mfe": 0.0,
        "bp_matrix": {},
        "per_base_unpaired": [0.5] * n,
        "n": n,
    }


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
    per_base_unpaired : list of float (1-indexed accessible probabilities from run_rnafold)
    window_start : int, 1-based CDS position of window start
    window_end : int, 1-based CDS position of window end
    cds_start_in_sequence : int
        0-based offset of the CDS start within the sequence passed to RNAfold.

    Returns
    -------
    float: mean accessibility in [0, 1]
    """
    seq_start = window_start + cds_start_in_sequence  # 1-indexed in sequence
    seq_end = window_end + cds_start_in_sequence

    values = []
    for pos in range(seq_start, seq_end + 1):
        if 1 <= pos <= len(per_base_unpaired):
            values.append(per_base_unpaired[pos - 1])  # per_base_unpaired is 0-indexed

    if not values:
        return 0.5
    return sum(values) / len(values)


if __name__ == "__main__":
    if check_rnafold_available():
        test_seq = "AAGACCGAGAGCAAGAAGGAGCGGCACGGCATGGCCATG"
        result = run_rnafold(test_seq)
        print(f"MFE structure: {result['mfe_structure']}")
        print(f"MFE: {result['mfe']} kcal/mol")
        print(f"Mean accessibility: {sum(result['per_base_unpaired'])/len(result['per_base_unpaired']):.3f}")
    else:
        print("RNAfold not available. Install ViennaRNA from https://www.tbi.univie.ac.at/RNA/")
