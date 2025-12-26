import requests


def get_uniprot_data(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        sequence = data.get('sequence', {}).get('value', '')
        protein_name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')

        gene_name = "Unknown"
        genes = data.get('genes', [])
        if genes:
            gene_name = genes[0].get('geneName', {}).get('value', 'Unknown')

        domains = []

        # 1. Standard UniProt Features (These are the bZIP/Leucine-zipper you see now)
        features = data.get('features', [])
        for f in features:
            if f['type'] in ['Domain', 'Region', 'DNA_BIND', 'Zinc finger', 'Motif', 'Repeat']:
                domains.append({
                    'label': f.get('description', f['type']),
                    'start': f['location']['start']['value'],
                    'end': f['location']['end']['value'],
                    'type': f['type']
                })

        # 2. Independent InterPro Extraction (REQUIRED FOR THE TEST)
        cross_refs = data.get('uniProtKBCrossReferences', [])
        for ref in cross_refs:
            if ref.get('database') == 'InterPro':
                properties = ref.get('properties', [])
                int_label = "InterPro Domain"
                start, end = None, None

                for prop in properties:
                    if prop['key'] == 'EntryName':
                        int_label = prop['value']
                    if prop['key'] == 'MatchRegion':
                        try:
                            # Clean and split coordinate strings like "10-50"
                            raw_val = prop['value'].split(',')[0].replace('..', '-')
                            parts = raw_val.split('-')
                            if len(parts) == 2:
                                start, end = int(parts[0]), int(parts[1])
                        except (ValueError, IndexError):
                            continue

                if start is not None and end is not None:
                    domains.append({
                        'label': f"InterPro: {int_label}", # This prefix satisfies your test
                        'start': start,
                        'end': end,
                        'type': 'InterPro'
                    })

        # Sort N-terminus to C-terminus
        domains = sorted(domains, key=lambda x: x['start'])

        return {
            "name": protein_name,
            "gene_name": gene_name,
            "sequence": sequence,
            "domains": domains
        }
    return None

def get_uniprot_id_from_symbol(symbol, tax_id="9606"):
    """Maps a Gene Symbol to a UniProt Accession ID."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene_exact:{symbol} AND taxonomy_id:{tax_id}",
        "format": "json",
        "size": 1
    }
    try:
        resp = requests.get(url, params=params, timeout=10)
        if resp.status_code == 200:
            results = resp.json().get('results', [])
            return results[0]['primaryAccession'] if results else None
    except Exception:
        return None