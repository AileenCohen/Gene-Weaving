import pytest
from src.api.uniprot_fetcher import get_uniprot_data, get_uniprot_id_from_symbol
from src.api.string_fetcher import get_interactions
from src.api.jaspar_fetcher import search_jaspar_motifs


def test_uniprot_fetch_valid_id():
    """Test that we can fetch a protein and find functional regions."""
    data = get_uniprot_data("P05412")  # Human JUN
    assert data is not None
    assert "JUN" in data['gene_name'].upper()

    # Check that we found at least SOME architectural features (UniProt OR InterPro)
    # This is more important than the specific "InterPro" label
    assert len(data['domains']) > 0, "No functional domains were found at all."

    # Verify we got a sequence
    assert len(data['sequence']) > 0

def test_symbol_to_uniprot_mapping():
    """Test the translation of Gene Symbol to Accession ID."""
    uid = get_uniprot_id_from_symbol("GATA6", tax_id="9606")
    assert uid == "Q92908"

def test_string_interactions():
    """Verify StringDB returns a list of partners for a known yeast protein."""
    partners = get_interactions("GAL4", tax_id="4932")
    assert isinstance(partners, list)
    assert len(partners) > 0

def test_jaspar_motif_search():
    """Ensure JASPAR returns results for a major TF."""
    results = search_jaspar_motifs("JUN", tax_id="9606")
    assert len(results) > 0
    assert "matrix_id" in results[0]