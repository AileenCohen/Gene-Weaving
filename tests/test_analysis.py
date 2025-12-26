import pytest
from src.analysis.primer_designer import get_optimized_dna, generate_primers
from src.analysis.crispr_designer import design_grnas


def test_codon_optimization():
    """Test that protein-to-DNA conversion maintains length (AA * 3)."""
    aa_seq = "METVAL"
    dna_human = get_optimized_dna(aa_seq, organism="Human")
    dna_yeast = get_optimized_dna(aa_seq, organism="Yeast")

    assert len(dna_human) == 18
    assert dna_human != dna_yeast  # Codon usage should differ


def test_primer_generation():
    """Verify primers contain the correct overhangs."""
    dna = "ATGCGTACGTAGCTAGCTAGCTAGCTAG"
    fwd_oh, rev_oh = "GAATTC", "GGATCC"
    primers = generate_primers(dna, fwd_oh, rev_oh)

    assert primers['fwd'].startswith(fwd_oh)
    assert primers['rev'].startswith(rev_oh)


def test_crispr_coordinate_logic():
    # Increase sequence length to 30 to give the 23bp window room to "breathe"
    test_dna = "G" * 30 + "AGG"
    guides = design_grnas(test_dna)
    assert len(guides) > 0
    assert "start_index" in guides[0]
    assert guides[0]['pam'] == "NGG"