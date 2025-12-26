

import requests
import pandas as pd


def search_jaspar_motifs(keyword: str, tax_id: str = "9606"):
    """
    Searches JASPAR matrices by name and taxonomy.
    """
    url = "https://jaspar.elixir.no/api/v1/matrix/"
    params = {
        "name": keyword,
        "tax_id": tax_id,
        "is_latest": True
    }

    try:
        resp = requests.get(url, params=params, timeout=10)
        if resp.status_code == 200:
            return resp.json().get('results', [])
        return []
    except Exception:
        return []


def get_pfm_data(matrix_id: str):
    """
    Fetches the Position Frequency Matrix (PFM) for a specific motif.
    """
    url = f"https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        pfm = resp.json().get('pfm')
        if pfm:
            return pd.DataFrame(pfm)
    except Exception:
        return None