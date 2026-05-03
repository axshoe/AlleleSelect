"""
accessibility.py
mRNA secondary structure accessibility scoring using RNAfold (Vienna RNA package).

RNAfold reference:
    Lorenz, R., et al. (2011). ViennaRNA Package 2.0. Algorithms Mol Biol 6(1):26.
    Free download: https://www.tbi.univie.ac.at/RNA/

Complete bug history:
  Bug 1: --noPS suppressed _dp.ps entirely. Removed.
  Bug 2: --id-prefix=target implies --auto-id which overrides the FASTA header and
         produces "target_0001_dp.ps", not "target_dp.ps". Removed that flag.
  Bug 3: Passing a file argument causes RNAfold on Windows to name output after the
         input *file* stem. Fixed by piping via stdin.

Definitive behaviour per RNAfold docs:
  - Plain text input (no FASTA header) -> output: dot.ps / rna.ps
  - FASTA input with ">myseq" header   -> output: myseq_dp.ps  (in cwd)
  - --id-prefix=X implies --auto-id, overrides FASTA header, produces X_0001_dp.ps

Correct approach: pipe ">target\n<seq>\n" via stdin, no --id-prefix, no --noPS.
Output will be target_dp.ps in work_dir. Full directory scan as fallback.
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
        RNA/DNA sequence (T will be converted to U).
    work_dir : str or None
        Directory where RNAfold writes output files. Uses tempdir if None.

    Returns
    -------
    dict with keys: mfe_structure, mfe, bp_matrix, per_base_unpaired, n
    """
    seq_rna = sequence.upper().replace("T", "U")
    n = len(seq_rna)

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="alleleselect_rnafold_")

    # Pipe FASTA via stdin. Key rules from official docs:
    # - FASTA header ">target" -> output file is "target_dp.ps" in cwd
    # - Do NOT use --id-prefix: it implies --auto-id which overrides the FASTA
    #   header and produces "target_0001_dp.ps" instead
    # - Do NOT use --noPS: suppresses all PostScript output including _dp.ps
    # - Piping via stdin (not file argument) avoids Windows path-stem naming
    fasta_input = f">target\n{seq_rna}\n"

    cmd = [
        "RNAfold",
        "-p",   # partition function + base-pair probability matrix (_dp.ps)
                # no --noPS, no --id-prefix
    ]

    try:
        result = subprocess.run(
            cmd,
            input=fasta_input,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=work_dir,   # RNAfold writes target_dp.ps here
        )
    except FileNotFoundError:
        warnings.warn(
            "RNAfold not found. Install ViennaRNA: https://www.tbi.univie.ac.at/RNA/\n"
            "Accessibility scores will be 0.5 (neutral) for all windows."
        )
        return _neutral_result(n)
    except subprocess.TimeoutExpired:
        warnings.warn("RNAfold timed out. Using neutral accessibility scores.")
        return _neutral_result(n)

    if result.returncode != 0:
        warnings.warn(
            f"RNAfold error (returncode {result.returncode}):\n{result.stderr[:500]}"
        )
        return _neutral_result(n)

    # Parse MFE structure and energy from stdout
    mfe_structure = ""
    mfe = 0.0
    for line in result.stdout.strip().split("\n"):
        m = re.match(r"^([.()\[\]{}|,]+)\s+\(\s*(-?\d+\.\d+)\)", line)
        if m:
            mfe_structure = m.group(1)
            mfe = float(m.group(2))
            break

    # Locate the dot-plot file.
    # Priority: target_dp.ps (expected), then scan directory for any *_dp.ps.
    # dot.ps is the fallback for old-style plain-text input (no FASTA header).
    candidate_dp_files = [
        os.path.join(work_dir, "target_dp.ps"),
        os.path.join(work_dir, "dot.ps"),
    ]
    try:
        for fname in sorted(os.listdir(work_dir)):
            full = os.path.join(work_dir, fname)
            if fname.endswith("_dp.ps") and full not in candidate_dp_files:
                candidate_dp_files.append(full)
    except OSError:
        pass

    bp_matrix = {}
    found_file = None
    for dp_file in candidate_dp_files:
        if os.path.exists(dp_file):
            bp_matrix = _parse_dp_file(dp_file)
            if bp_matrix:
                found_file = dp_file
                break

    if not bp_matrix:
        warnings.warn(
            "RNAfold completed but no usable _dp.ps was found.\n"
            "Accessibility scores will be 0.5 (neutral) for all windows.\n"
            f"  work_dir : {work_dir}\n"
            f"  files    : {os.listdir(work_dir)}\n"
            "If a *_dp.ps file appears above, open it and check it contains\n"
            "lines of the form: <int> <int> <float> ubox"
        )
        return {
            "mfe_structure": mfe_structure,
            "mfe": mfe,
            "bp_matrix": {},
            "per_base_unpaired": [0.5] * n,
            "n": n,
        }

    # Compute per-base unpaired probability.
    # RNAfold dp.ps stores sqrt(P(i,j)) in ubox entries; square to get true prob.
    # P_unpaired(i) = 1 - sum_j P(i,j)
    per_base_paired = [0.0] * (n + 1)   # 1-indexed; [0] unused
    for (i, j), sqrt_prob in bp_matrix.items():
        prob = sqrt_prob ** 2
        if 1 <= i <= n:
            per_base_paired[i] += prob
        if 1 <= j <= n:
            per_base_paired[j] += prob

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
    Parse RNAfold _dp.ps PostScript file.
    Returns {(i, j): sqrt_probability} for all ubox entries where i < j.
    """
    bp_matrix = {}
    pattern = re.compile(r"(\d+)\s+(\d+)\s+([\d.e+\-]+)\s+ubox")
    try:
        with open(dp_path, "r", errors="replace") as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    i, j, sq = int(m.group(1)), int(m.group(2)), float(m.group(3))
                    if i < j:
                        bp_matrix[(i, j)] = sq
    except Exception as e:
        warnings.warn(f"Could not parse {dp_path}: {e}")
    return bp_matrix


def _neutral_result(n: int) -> dict:
    """Fallback when RNAfold is unavailable: all positions 0.5 (neutral)."""
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
    Mean unpaired probability across an ASO candidate window.

    Parameters
    ----------
    per_base_unpaired : 0-indexed list from run_rnafold
    window_start : 1-based CDS position of window start
    window_end   : 1-based CDS position of window end
    cds_start_in_sequence : 0-based offset of CDS start inside the RNAfold sequence
        e.g. if RNAfold ran on a ±200 nt window starting at CDS pos 375,
        this is 374 (= 375 - 1).
    """
    seq_start = window_start + cds_start_in_sequence
    seq_end   = window_end   + cds_start_in_sequence
    values = [
        per_base_unpaired[pos - 1]
        for pos in range(seq_start, seq_end + 1)
        if 0 <= pos - 1 < len(per_base_unpaired)
    ]
    return sum(values) / len(values) if values else 0.5


if __name__ == "__main__":
    if check_rnafold_available():
        test_seq = "AAGACCGAGAGCAAGAAGGAGCGGCACGGCATGGCCATG"
        wd = tempfile.mkdtemp(prefix="alleleselect_test_")
        r = run_rnafold(test_seq, work_dir=wd)
        print(f"Length      : {r['n']}")
        print(f"MFE struct  : {r['mfe_structure']}")
        print(f"MFE energy  : {r['mfe']} kcal/mol")
        print(f"BP pairs    : {len(r['bp_matrix'])}")
        mean_acc = sum(r['per_base_unpaired']) / len(r['per_base_unpaired'])
        print(f"Mean access : {mean_acc:.3f}")
        if not r["bp_matrix"]:
            print(f"WARNING: bp_matrix is empty. Check {wd} for files.")
        else:
            print("Accessibility scoring is working correctly.")
    else:
        print("RNAfold not found. Install from https://www.tbi.univie.ac.at/RNA/")