"""
hgvs_parser.py
Parses HGVS coding DNA notation (e.g. c.575G>A) and validates against a CDS sequence.
"""

import re


class HGVSParseError(Exception):
    pass


def parse_hgvs_coding(hgvs_str: str) -> dict:
    """
    Parse a coding-DNA HGVS substitution string.

    Parameters
    ----------
    hgvs_str : str
        e.g. "c.575G>A"

    Returns
    -------
    dict with keys: position (int, 1-based), ref (str), alt (str)

    Raises
    ------
    HGVSParseError on malformed input.
    """
    pattern = r"^c\.(\d+)([ACGTacgt])>([ACGTacgt])$"
    m = re.match(pattern, hgvs_str.strip())
    if not m:
        raise HGVSParseError(
            f"Cannot parse HGVS string '{hgvs_str}'. "
            f"Expected format: c.<position><ref>><alt> (e.g. c.575G>A)"
        )
    position = int(m.group(1))
    ref = m.group(2).upper()
    alt = m.group(3).upper()
    if ref == alt:
        raise HGVSParseError(f"ref and alt are identical: '{ref}'")
    return {"position": position, "ref": ref, "alt": alt}


def validate_against_cds(parsed: dict, cds_sequence: str) -> None:
    """
    Confirm that the reference nucleotide at the stated position matches the CDS.

    Parameters
    ----------
    parsed : dict returned by parse_hgvs_coding
    cds_sequence : str, the full CDS (coding) sequence (1-indexed in biology, 0-indexed in Python)

    Raises
    ------
    HGVSParseError if there is a mismatch.
    """
    pos = parsed["position"]
    if pos < 1 or pos > len(cds_sequence):
        raise HGVSParseError(
            f"Position {pos} is out of range for CDS of length {len(cds_sequence)}."
        )
    actual_base = cds_sequence[pos - 1].upper()
    if actual_base != parsed["ref"]:
        raise HGVSParseError(
            f"CDS base at position {pos} is '{actual_base}', "
            f"but HGVS ref says '{parsed['ref']}'. Check transcript ID or HGVS notation."
        )


def apply_variant_to_sequence(cds_sequence: str, parsed: dict) -> str:
    """
    Return a mutant CDS with the alt nucleotide substituted at the given position.
    """
    pos = parsed["position"]
    mutant = (
        cds_sequence[: pos - 1] + parsed["alt"] + cds_sequence[pos:]
    )
    return mutant


if __name__ == "__main__":
    # Quick demo for R192Q
    test_hgvs = "c.575G>A"
    result = parse_hgvs_coding(test_hgvs)
    print(f"Parsed: position={result['position']}, ref={result['ref']}, alt={result['alt']}")
