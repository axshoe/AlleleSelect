"""
fetcher.py
Retrieves CACNA1A mRNA/CDS sequence from Ensembl REST API by transcript ID.
Also generates mutant sequence given a parsed HGVS dict.
"""

import requests
import time
from alleleselect.sequence.hgvs_parser import apply_variant_to_sequence

ENSEMBL_REST_BASE = "https://rest.ensembl.org"
DEFAULT_TRANSCRIPT = "ENST00000360228"  # CACNA1A canonical transcript


class EnsemblFetchError(Exception):
    pass


def fetch_cds_sequence(transcript_id: str = DEFAULT_TRANSCRIPT, retries: int = 3) -> str:
    """
    Fetch the CDS (coding sequence) for a given Ensembl transcript ID.

    Parameters
    ----------
    transcript_id : str
        Ensembl transcript ID, with or without version suffix (e.g. ENST00000360228 or ENST00000360228.10)
    retries : int
        Number of retry attempts on network failure.

    Returns
    -------
    str : CDS nucleotide sequence (A/C/G/T, uppercase)
    """
    # Strip version suffix for the API call
    base_id = transcript_id.split(".")[0]
    url = f"{ENSEMBL_REST_BASE}/sequence/id/{base_id}"
    params = {"content-type": "application/json", "type": "cds"}
    headers = {"Content-Type": "application/json"}

    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                seq = data.get("seq", "").upper()
                if not seq:
                    raise EnsemblFetchError("Empty sequence returned from Ensembl.")
                return seq
            elif resp.status_code == 429:
                wait = int(resp.headers.get("Retry-After", 5))
                print(f"Rate limited. Waiting {wait}s...")
                time.sleep(wait)
            else:
                raise EnsemblFetchError(
                    f"Ensembl returned status {resp.status_code}: {resp.text[:200]}"
                )
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                raise EnsemblFetchError(f"Network error fetching sequence: {e}")

    raise EnsemblFetchError("Max retries exceeded fetching CDS sequence.")


def fetch_mrna_sequence(transcript_id: str = DEFAULT_TRANSCRIPT) -> str:
    """
    Fetch the full mRNA sequence (includes UTRs) for secondary structure analysis.
    """
    base_id = transcript_id.split(".")[0]
    url = f"{ENSEMBL_REST_BASE}/sequence/id/{base_id}"
    params = {"content-type": "application/json", "type": "cdna"}
    resp = requests.get(url, params=params, timeout=30)
    if resp.status_code != 200:
        raise EnsemblFetchError(f"mRNA fetch failed: {resp.status_code}")
    return resp.json().get("seq", "").upper()


def generate_mutant_cds(wildtype_cds: str, parsed_hgvs: dict) -> str:
    """
    Apply parsed HGVS variant to the wildtype CDS to produce mutant CDS.
    """
    return apply_variant_to_sequence(wildtype_cds, parsed_hgvs)


def extract_window(sequence: str, center_pos: int, flank: int = 200) -> tuple:
    """
    Extract a sequence window of ±flank nucleotides around a center position (1-based).
    Returns (window_sequence, window_start_1based).
    """
    pos0 = center_pos - 1  # convert to 0-indexed
    start = max(0, pos0 - flank)
    end = min(len(sequence), pos0 + flank + 1)
    return sequence[start:end], start + 1  # return 1-based start


if __name__ == "__main__":
    print("Fetching CACNA1A CDS from Ensembl...")
    cds = fetch_cds_sequence()
    print(f"CDS length: {len(cds)} nt")
    print(f"Bases 570-580 (around R192Q site): {cds[569:580]}")
