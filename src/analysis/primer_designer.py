def get_optimized_dna(protein_seq, organism="Human"):
    # High-efficiency codon tables
    tables = {
        "Human": {'A':'GCC','C':'TGC','D':'GAC','E':'GAG','F':'TTC','G':'GGC','H':'CAC','I':'ATC','K':'AAG','L':'CTG','M':'ATG','N':'AAC','P':'CCC','Q':'CAG','R':'CGC','S':'TCC','T':'ACC','V':'GTG','W':'TGG','Y':'TAC','*':'TAA'},
        "Yeast": {'A':'GCT','C':'TGT','D':'GAT','E':'GAA','F':'TTT','G':'GGT','H':'CAT','I':'ATT','K':'AAA','L':'TTG','M':'ATG','N':'AAT','P':'CCA','Q':'CAA','R':'AGA','S':'TCT','T':'ACT','V':'GTT','W':'TGG','Y':'TAT','*':'TAA'}
    }
    table = tables.get(organism, tables["Human"])
    return "".join([table.get(aa, 'NNN') for aa in protein_seq])

def scan_restriction_sites(dna):
    enzymes = {"EcoRI":"GAATTC", "BamHI":"GGATCC", "HindIII":"AAGCTT", "NotI":"GCGGCCGC", "XhoI":"CTCGAG"}
    return [name for name, site in enzymes.items() if site in dna.upper()]

def calculate_tm(seq):
    return (2 * (seq.count('A') + seq.count('T'))) + (4 * (seq.count('G') + seq.count('C')))

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def generate_primers(dna, ov_f, ov_r):
    bf, br = dna[:20], reverse_complement(dna[-20:])
    return {"fwd": f"{ov_f}{bf}", "rev": f"{ov_r}{br}", "fwd_bind": bf, "rev_bind": br}

def get_primer_quality(full_seq, binding_part):
    tm = calculate_tm(binding_part)
    gc = (full_seq.count('G') + full_seq.count('C')) / len(full_seq) * 100
    status = "✅" if 40 <= gc <= 60 else "⚠️"
    return {"tm": tm, "gc": gc, "status": status, "warnings": ""}