import streamlit as st
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logomaker
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Path setup
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from api.uniprot_fetcher import get_uniprot_data, get_uniprot_id_from_symbol
from analysis.disorder_analyzer import get_disorder_scores, analyze_disorder
from api.jaspar_fetcher import search_jaspar_motifs, get_pfm_data
from analysis.primer_designer import (
    get_optimized_dna,
    generate_primers,
    get_primer_quality,
    scan_restriction_sites
)
from analysis.crispr_designer import design_grnas
from api.string_fetcher import get_interactions


# --- VISUALIZATION: PROTEIN ARCHITECTURE ---
def plot_protein_architecture(total_len, domains, idrs):
    """Draws a horizontal diagram of domains and IDRs with full labels."""
    fig, ax = plt.subplots(figsize=(12, 3.0))  # Increased height for labels
    # Backbone
    ax.add_patch(patches.Rectangle((1, 0.8), total_len, 0.4, color='lightgray', alpha=0.3))

    # Plot IDRs (Thin orange)
    for i in idrs:
        ax.add_patch(patches.Rectangle((i['start'], 0.95), i['end'] - i['start'], 0.1, color='orange', alpha=0.6))

    # Plot Domains (Colored blocks)
    colors = ['#4A90E2', '#50E3C2', '#F5A623', '#D0021B', '#8B572A']
    for idx, d in enumerate(domains):
        color = '#9B51E0' if "InterPro" in d.get('label', '') else colors[idx % len(colors)]
        ax.add_patch(patches.Rectangle((d['start'], 0.7), d['end'] - d['start'], 0.6,
                                       edgecolor='black', facecolor=color, alpha=0.8))

        # Add labels to larger domains
        if (d['end'] - d['start']) > total_len * 0.02:
            label_text = d['label'].split(':')[-1].strip()

            ax.text(
                (d['start'] + d['end']) / 2,
                1.5,
                label_text,
                ha='left',
                va='bottom',
                fontsize=9,
                fontweight='bold',
                rotation=35  # Better angle for long text
            )

    ax.set_xlim(0, total_len + 10)
    ax.set_ylim(0, 4)  # Generous headroom for labels
    ax.axis('off')
    st.pyplot(fig)


# --- CACHE ---
@st.cache_data(show_spinner=False)
def get_cached_analysis(sequence):
    return get_disorder_scores(sequence), analyze_disorder(sequence)


# --- APP SETUP ---
st.set_page_config(page_title="Gene Weaving", layout="wide", page_icon="🧬")
st.title("🧬 Gene Weaving: TF Designer")

# Session State
if "history" not in st.session_state: st.session_state.history = []
if "protein_data" not in st.session_state: st.session_state.protein_data = None
if "primer_list" not in st.session_state: st.session_state.primer_list = []
if "start_val" not in st.session_state: st.session_state.start_val = 1
if "end_val" not in st.session_state: st.session_state.end_val = 100

with st.sidebar:
    st.header("Search Settings")
    u_id = st.text_input("UniProt ID:", value=st.session_state.get('current_id', "P01106")).strip().upper()
    organism = st.radio("Target Host:", ["Human", "Yeast"])
    tax_id = "9606" if organism == "Human" else "4932"

    if st.button("Analyze Protein"):
        data = get_uniprot_data(u_id)
        if data:
            st.session_state.protein_data = data
            st.session_state.current_id = u_id
            st.session_state.start_val, st.session_state.end_val = 1, len(data['sequence'])
            st.session_state.motifs = search_jaspar_motifs(data.get('gene_name', 'Unknown'), tax_id=tax_id)
            if u_id not in st.session_state.history: st.session_state.history.append(u_id)
            st.rerun()
        else:
            st.error("Protein not found.")

    if st.session_state.history:
        st.write("🕒 History")
        for h_id in st.session_state.history[-3:]:
            if st.button(f"Go to {h_id}", key=f"hist_nav_{h_id}"):
                st.session_state.current_id = h_id
                st.session_state.protein_data = get_uniprot_data(h_id)
                st.rerun()

# --- MAIN INTERFACE ---
if st.session_state.protein_data:
    data = st.session_state.protein_data
    scores, idrs = get_cached_analysis(data['sequence'])

    # Header
    st.header(f"Protein: {data['name']} ({st.session_state.current_id})")

    # Architecture Visual
    st.subheader("📊 Architecture & Disorder")

    plot_protein_architecture(len(data['sequence']), data['domains'], idrs)
    st.line_chart(pd.DataFrame({"Disorder Probability": scores}), color="#FF4B4B")

    tab_design, tab_motifs, tab_interact, tab_crispr = st.tabs([
        "🏗️ Construct Designer", "🧬 Binding Motifs", "🤝 Interactions", "✂️ CRISPR/gRNA"
    ])

    with tab_design:
        col_view, col_input = st.columns([1, 1])
        with col_view:
            st.subheader("📋 Structural Features")
            if data.get('domains'):
                for d in data['domains']:
                    c_info, c_btn = st.columns([4, 1])
                    c_info.info(f"**{d['type']}**: {d['label']} ({d['start']}-{d['end']})")
                    # UNIQUE KEY FIX
                    btn_key = f"set_{d['type']}_{d['start']}_{d['end']}"
                    if c_btn.button("Set", key=btn_key):
                        st.session_state.start_val, st.session_state.end_val = d['start'], d['end']
                        st.rerun()

        with col_input:
            st.subheader("🛠️ Parameters")
            cur_s = st.number_input("Start Residue:", 1, len(data['sequence']), value=int(st.session_state.start_val))
            cur_e = st.number_input("End Residue:", 1, len(data['sequence']), value=int(st.session_state.end_val))

            # Sync session state
            st.session_state.start_val, st.session_state.end_val = cur_s, cur_e

            target_aa = data['sequence'][int(cur_s) - 1: int(cur_e)]
            target_dna = get_optimized_dna(target_aa, organism=organism)

            if target_aa:
                pa = ProteinAnalysis(target_aa)
                st.write(f"**MW:** {pa.molecular_weight() / 1000:.1f} kDa | **pI:** {pa.isoelectric_point():.1f}")

            sites = scan_restriction_sites(target_dna)
            if sites:
                st.error(f"⚠️ Internal Restriction Sites: {', '.join(sites)}")
            else:
                st.success("✅ No internal restriction sites.")

            if st.button("Add Construct to Collection"):
                p = generate_primers(target_dna, "GAATTC", "GGATCC")
                st.session_state.primer_list.append(
                    {"Label": f"{st.session_state.current_id}_{cur_s}-{cur_e}", "Fwd": p['fwd'], "Rev": p['rev']})
                st.toast("Construct Saved!")

    with tab_motifs:
        if st.session_state.get('motifs'):
            m_names = [m['name'] for m in st.session_state.motifs]
            sel = st.selectbox("Select Motif Version:", m_names)
            m_id = st.session_state.motifs[m_names.index(sel)]['matrix_id']
            pfm = get_pfm_data(m_id)
            if pfm is not None:
                fig_m, ax_m = plt.subplots(figsize=(8, 2.5))
                logomaker.Logo(logomaker.transform_matrix(pfm, from_type='counts', to_type='information'), ax=ax_m)
                st.pyplot(fig_m)
        else:
            st.warning("No direct DNA-binding motifs found.")
            st.info("💡 Hint: This may be a co-activator. Check 'Interactions' for partners.")

    with tab_interact:
        st.subheader("🤝 Physical Interaction Network (String-DB)")

        partners = get_interactions(data.get('gene_name', 'Unknown'), tax_id=tax_id)
        if partners:
            for p in partners:
                c1, c2 = st.columns([3, 1])
                c1.code(p)
                if c2.button("Analyze Partner", key=f"switch_{p}"):
                    new_id = get_uniprot_id_from_symbol(p, tax_id=tax_id)
                    if new_id:
                        st.session_state.current_id = new_id
                        st.session_state.protein_data = None
                        st.rerun()
        else:
            st.write("No high-confidence partners found.")

    with tab_crispr:
        st.subheader("✂️ gRNA Domain Mapping")

        if st.button("Run CRISPR Scan"):
            grnas = design_grnas(target_dna)
            if grnas:
                mapping = []
                for g in grnas:
                    abs_pos = cur_s + (g['start_index'] // 3)
                    loc = "Linker Region"
                    for d in data['domains']:
                        if d['start'] <= abs_pos <= d['end']:
                            loc = f"🎯 Hits {d['label']}"
                    mapping.append({"gRNA ID": g['label'], "Abs Pos (AA)": int(abs_pos), "Location": loc,
                                    "GC%": f"{g['gc']:.0f}%"})
                st.table(pd.DataFrame(mapping))
            else:
                st.error("No NGG PAM sites found in this selection.")

    if st.session_state.primer_list:
        st.divider()
        st.subheader("📋 Saved Constructs")
        st.dataframe(pd.DataFrame(st.session_state.primer_list), use_container_width=True)
        csv = pd.DataFrame(st.session_state.primer_list).to_csv(index=False).encode('utf-8')
        st.download_button("📥 Download Primer CSV", csv, f"primers_{st.session_state.current_id}.csv", "text/csv")
else:
    st.info("Enter a UniProt ID to begin (e.g., P01106 for Human MYC).")