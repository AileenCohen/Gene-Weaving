import requests


def get_interactions(gene_name, tax_id="9606"):
    """
    Fetches protein-protein interactions from StringDB.
    9606 = Human, 4932 = Yeast
    """
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": gene_name,
        "species": tax_id,
        "limit": 5
    }
    try:
        resp = requests.get(url, params=params, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            partners = set()
            for entry in data:
                # We want the name of the partner protein
                partners.add(entry['preferredName_B'])

            # Filter out the original protein from the list
            if gene_name.upper() in [p.upper() for p in partners]:
                partners = {p for p in partners if p.upper() != gene_name.upper()}

            return list(partners)
    except Exception:
        return []