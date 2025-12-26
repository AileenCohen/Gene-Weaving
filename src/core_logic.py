import sys
import os

# This line tells Python to look in the current folder for the other files
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Now we can import your tools
from api.uniprot_fetcher import get_uniprot_data
from analysis.disorder_analyzer import analyze_disorder, get_disorder_scores


def weave_protein_report(uniprot_id):
    # 1. Get data from UniProt
    data = get_uniprot_data(uniprot_id)
    if not data:
        return

    print(f"\n--- Gene Weaving Report for {data['name']} ---")

    # 2. Run Disorder Analysis
    idrs = analyze_disorder(data['sequence'])

    # 3. Combine and Print a Summary
    print("\n[Architecture Summary]")

    # Show structural domains from UniProt
    print("\nKnown Structural Domains:")
    for d in data['domains']:
        # We check for 'Domain' or 'Region' since TFs use both labels
        if d['type'] in ['Domain', 'Region']:
            print(f"  -> {d['label']}: residues {d['start']}-{d['end']}")

    # Show predicted floppy regions
    print("\nPredicted Disordered Regions (IDRs):")
    if not idrs:
        print("  (No large IDRs met the automatic threshold)")
    else:
        for idr in idrs:
            print(f"  ~ IDR: residues {idr['start']}-{idr['end']}")

    # Optional: Print the average disorder score to see if it's "mostly floppy"
    scores = get_disorder_scores(data['sequence'])
    avg_disorder = sum(scores) / len(scores)
    print(f"\nAverage Protein Disorder Score: {avg_disorder:.2f}")


if __name__ == "__main__":
    # P01106 = MYC (Very disordered)
    # P10275 = Androgen Receptor (Half disordered)
    weave_protein_report("P01106")