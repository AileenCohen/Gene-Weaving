---
layout: default
---

<link rel="stylesheet" href="{{ '/assets/css/main.css' | relative_url }}">

<p align="center">
  <img src="assets/img/Gene_Weaving.jpg" width="300" alt="Gene Weaving Logo">
</p>

# Gene Weaving: TF Designer

**A specialized bioinformatics workbench for Transcription Factor engineering and construct design.**

## Overview
**Gene Weaving** is a lightweight Python toolkit and interactive workbench designed to streamline the design of transcription factor (TF) constructs. It bridges the gap between raw sequence data and lab-ready cloning strategies by centralizing protein architecture, disorder prediction, and primer design into a single automated workflow.

The project focuses on the unique challenges of **Intrinsically Disordered Proteins (IDPs)** and complex TF assemblies, providing researchers with a visual and analytical "loom" to weave together domains, disordered regions, and binding sites.

---

## Motivation
Designing constructs for TFs often requires manual data collection from multiple siloed sources. Gene Weaving aims to:
* **Centralize Annotations:** Aggregate data from UniProt, JASPAR, and String-DB.
* **Visualize "Hidden" Architecture:** Map InterPro structural domains alongside **MetaPredict** disorder profiles to identify surgical boundaries for truncations.
* **Automate Interaction Discovery:** Identify recruitment partners for co-activators that lack direct DNA-binding motifs.
* **Simplify Precision Engineering:** Automate CRISPR/gRNA mapping and primer generation for deletions or substitutions.

---

## Features

### 1. Integrated Data Retrieval
* **UniProt & InterPro Mapping:** Fetches primary sequences and maps manually curated features alongside signature-based structural domains.
* **JASPAR Binding Motifs:** Retrieves and renders DNA-binding sequence logos.
* **String-DB Interaction Network:** Identifies high-confidence physical interaction partners with a "one-click" analysis feature.

### 2. Disorder & Architecture Visualization
* **MetaPredict Integration:** Provides residue-level disorder probability scores to identify Intrinsically Disordered Regions (IDRs).
* **Protein Mapping:** Generates a to-scale diagram of the protein, overlaying IDRs and folded domains.

### 3. Precision Design Tools
* **Automated Primer Designer:** Generates sequences with customizable 5'/3' overhangs (e.g., for Gibson Assembly).
* **CRISPR/gRNA Mapper:** Scans regions for NGG PAM sites and maps gRNAs back to annotated domains to ensure functional targeting.
* **Cloning Validation:** Screens synthetic DNA for internal restriction sites to prevent experimental failures.

---

## Project Structure
```text
GeneWeaving/
├── src/
│   ├── api/        # Interfaces for UniProt, JASPAR, and String-DB
│   ├── analysis/   # Disorder (MetaPredict), CRISPR (gRNA), and Primer logic
│   └── gui/        # Streamlit-based interactive workbench
├── tests/          # Pytest suite for API and analysis validation
├── requirements.txt
└── .gitignore      # Optimized for Python and Bioinformatics tools
```
---

## Installation & Setup
1. Clone the repository:
   ```
   git clone https://github.com/AileenCohen/Gene-Weaving.git
   cd Gene-Weaving
   ```
2. Create a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```
4. Run the app:
    ```
    python main.py
    ```
## Testing:
To ensure coordinate math and API connections are functioning correctly, run the automated test suite:
```
python -m pytest
```
---
## Example Workflow
- Search: Input a UniProt ID (e.g., P05412 for Human JUN).
- Explore: Visualize the bZIP domain and flanking IDRs on the Protein Map.
- Interact: View partners like FOS or ATF2 and switch to them instantly.
- Design: Highlight a domain to generate gRNAs for a knockout or primers for a truncation.
- Export: Download a CSV of primers and coordinates for bench use.
---
## Dependencies
- Python 3.10+
- Streamlit (UI)
- Biopython (Sequence handling)
- MetaPredict (IDR prediction)
- Logomaker (Motif visualization)
- Requests (API management)
- Pytest (Validation)
  
---

Gene Weaving is a project focused on providing accessible, practical tools for transcription factor construct design.  
Its scope is intentionally modest, emphasizing usability and core functionality over full automation.
If you have any additional ideas, please open an issue and let me know.

Created by Aileen Cohen.
