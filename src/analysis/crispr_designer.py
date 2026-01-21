def design_grnas(dna_seq):
    """
    Finds NGG PAM sites and returns gRNA candidates with coordinates.
    """
    results = []
    # Scan DNA for PAM (NGG)
    for i in range(len(dna_seq) - 23):
        # Check if the 21st and 22nd bases are 'GG'
        if dna_seq[i + 21:i + 23] == "GG":
            guide = dna_seq[i:i + 20]
            gc = (guide.count('G') + guide.count('C')) / 20 * 100

            results.append({
                "label": f"gRNA_{i + 1}",
                "sequence": guide,
                "pam": "NGG",
                "gc": gc,
                "start_index": i  
            })

    return results[:5]
