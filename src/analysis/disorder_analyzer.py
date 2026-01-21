import metapredict as mp


def get_disorder_scores(sequence: str):
    """Returns the per-residue disorder scores using metapredict."""
    if not sequence:
        return []
    return mp.predict_disorder(sequence)


def analyze_disorder(sequence: str):
    """
    Finds IDRs by checking metapredict scores against a threshold.
    Returns 1-based coordinates for UniProt compatibility.
    """
    if not sequence:
        return []

    scores = mp.predict_disorder(sequence)
    threshold = 0.5
    min_length = 30

    idrs = []
    current_idr_start = None

    for i, score in enumerate(scores):
        if score >= threshold:
            if current_idr_start is None:
                current_idr_start = i + 1  
        else:
            if current_idr_start is not None:
                current_idr_end = i
                if (current_idr_end - current_idr_start + 1) >= min_length:
                    idrs.append({
                        "start": current_idr_start,
                        "end": current_idr_end,
                        "type": "IDR"
                    })
                current_idr_start = None

    # Handle C-terminal IDR
    if current_idr_start is not None:
        current_idr_end = len(scores)
        if (current_idr_end - current_idr_start + 1) >= min_length:
            idrs.append({
                "start": current_idr_start,
                "end": current_idr_end,
                "type": "IDR"
            })


    return idrs
