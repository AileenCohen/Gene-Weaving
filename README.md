# Gene Weaving  
A lightweight Python toolkit and GUI concept for assisting with transcription factor construct design

## Overview
**Gene Weaving** is a developing bioinformatics tool aimed at simplifying common sequence-level design tasks for transcription factor (TF) research.  
Its goal is to provide a single interface where researchers can fetch basic TF annotations, preview domain/disorder architecture, and generate simple primer designs for common cloning workflows.

The project is intended as a practical assistant for small-scale TF engineering tasks, with an emphasis on clarity, reproducibility and automation.

---

## Motivation
Designing constructs for transcription factors often requires collecting information from several sources: UniProt for sequences, InterPro for domain boundaries, disorder prediction tools for IDR regions, and then manually incorporating that information into cloning strategies.  
Although none of these steps are individually complex, together they can slow down experimental planning.

Gene Weaver aims to:
- centralize essential annotation steps,
- offer straightforward primer design for basic constructs,
- and provide organized outputs suitable for lab documentation.

The goal is not to replace full-featured bioinformatics platforms, but to streamline routine design work for TF-focused studies.

---

## Project Structure (Draft)
```
GeneWeaver/
│
├─ src/
│ ├─ api/ # Basic interfaces for external data sources
│ ├─ analysis/ # Sequence parsing and annotation utilities
│ ├─ primer_design/ # Simple primer design functions
│ └─ gui/ # Early-stage GUI prototype (Streamlit or PyQt)
│
├─ examples/
├─ tests/
├─ README.md
└─ requirements.txt
```

---

## Features (Current & Planned Scope)

### Data Retrieval (Initial Functionality)
- Retrieve protein sequences from UniProt
- Fetch binding motifs (if available) from JASPAR  
- Load domain annotations from InterPro  
- Run basic disorder predictions (MetaPredict integration planned)  

### Basic Construct & Primer Design
- Generate primers for:
  - small deletions  
  - N-/C-terminal truncations  
  - simple substitution designs  
- Provide standard overhang options suitable for Gibson Assembly  

These design tools focus on small, clearly defined edits rather than complex multi-fragment workflows.

### Validation (Limited)
- Estimate melting temperature (Tm)  
- Calculate GC%  
- Screen for obvious hairpins or dimers  

### Export
- Primer lists to CSV  
- Basic FASTA export for modified sequences  

---

## Example Workflow
1. User inputs a UniProt ID (e.g., `P10275` for human MYC).  
2. Gene Weaver retrieves:
   - primary sequence  
   - InterPro-annotated domain boundaries  
   - predicted disorder regions (if enabled)  
3. The user selects a simple design (e.g., remove an annotated IDR or create a small truncation).  
4. The toolkit generates basic primers and provides Tm/GC estimates.  
5. Outputs are exported to CSV or FASTA for bench use.

This workflow is intended for straightforward construct design rather than comprehensive mutagenesis pipelines.

---

## GUI (Prototype Status)
A minimal GUI is planned to:
- accept UniProt input  
- display annotated sequence regions  
- offer primer design options  
- export results  

The interface will prioritize clarity and ease of use over feature complexity.

---

## Dependencies
- **Python 3.8+**  
- `requests` — external data retrieval  
- `Biopython` — sequence utilities  
- `pandas` — CSV handling  
- `metapredict` — disorder prediction  
- `Streamlit` or `PyQt` — basic GUI  
- `pytest` — testing  

---

Gene Weaver is a project focused on providing accessible, practical tools for transcription factor construct design.  
Its scope is intentionally modest, emphasizing usability and core functionality over full automation.
If you have any additional ideas, please open an issue and let me know.
