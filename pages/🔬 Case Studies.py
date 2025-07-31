import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------- UI Styling -----------------
st.set_page_config(page_title="Biological Case Studies", layout="wide")

st.markdown("""
<style>
body {
    background-color: #f7fafc;
    font-family: 'Segoe UI', sans-serif;
    color: #333;
}
h1, h2, h3 {
    color: #1e3a8a;
}
.stTabs [role="tablist"] {
    background-color: #e0f2fe;
    border-radius: 0.5rem;
    padding: 0.5rem
}
.stTabs [role="tab"][aria-selected="true"] {
    background-color: #1e40af !important;
    color: white !important;
    border-radius: 0.5rem;
    padding: 0.5rem
}
.stButton>button {
    background-color: #1d4ed8;
    color: white;
    font-weight: bold;
    border-radius: 6px;
    padding: 0.5rem 1rem;
}
.stButton>button:hover {
    background-color: #1e3a8a;
    transform: scale(1.02);
}
</style>
""", unsafe_allow_html=True)

# ----------------- Constants -----------------
MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_PENALTY = -2

# ----------------- Alignment Functions -----------------
def needleman_wunsch(ref, pat, match=MATCH_SCORE, mismatch=MISMATCH_SCORE, gap=GAP_PENALTY):
    m, n = len(ref), len(pat)
    score = np.zeros((m+1, n+1), dtype=int)
    for i in range(1, m+1): score[i][0] = i * gap
    for j in range(1, n+1): score[0][j] = j * gap
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = score[i-1][j-1] + (match if ref[i-1] == pat[j-1] else mismatch)
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap
            score[i][j] = max(diag, delete, insert)
    return score

def traceback(ref, pat, S, match=MATCH_SCORE, mismatch=MISMATCH_SCORE, gap=GAP_PENALTY):
    i, j = len(ref), len(pat)
    aligned_ref, aligned_pat = "", ""
    while i > 0 and j > 0:
        score = S[i][j]
        if score == S[i-1][j-1] + (match if ref[i-1] == pat[j-1] else mismatch):
            aligned_ref = ref[i-1] + aligned_ref
            aligned_pat = pat[j-1] + aligned_pat
            i -= 1; j -= 1
        elif score == S[i-1][j] + gap:
            aligned_ref = ref[i-1] + aligned_ref
            aligned_pat = '-' + aligned_pat
            i -= 1
        else:
            aligned_ref = '-' + aligned_ref
            aligned_pat = pat[j-1] + aligned_pat
            j -= 1
    while i > 0:
        aligned_ref = ref[i-1] + aligned_ref
        aligned_pat = '-' + aligned_pat
        i -= 1
    while j > 0:
        aligned_ref = '-' + aligned_ref
        aligned_pat = pat[j-1] + aligned_pat
        j -= 1
    return aligned_ref, aligned_pat

def highlight_inline_differences(ref_seq, mut_seq):
    result = ""
    for r, m in zip(ref_seq, mut_seq):
        if r == m:
            result += m
        elif r == "-":
            result += f"<span style='color: blue'><b>{m}</b></span>"
        elif m == "-":
            result += f"<span style='color: gray'><b>{r}</b></span>"
        else:
            result += f"<span style='color: red'><b>{m}</b></span>"
    return f"<code style='font-size: 18px; font-family: monospace'>{result}</code>"

# ------------------ Title & Tabs ------------------
st.title("🔬 Case Studies")
mutation_tabs = st.tabs(["🧬 Evolution (Cytochrome C)", "🧪 Mutation Case Studies"])

# ------------------ Tab 1: Evolution ------------------
with mutation_tabs[0]:
    st.header("🌍 Evolution of Primates via Cytochrome C Protein")
    st.markdown("Cytochrome C is highly conserved across species, ideal for studying evolutionary divergence.")

    cytochrome_data = {
        "Species": ["Human", "Chimpanzee", "Gorilla", "Orangutan", "Macaque"],
        "Sequence": [
            "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKEERADLIAYLKE",
            "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIYAGIKKKEERADLIAYLKE",
            "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIYAGIKKKEERADLIAYLKE",
            "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKAGQAPGFSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIYAGIKKKEERADLIAYLKE",
            "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKAGQAPGFSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIYAGIKKKEERADLIAYLKE"
        ]
    }

    df = pd.DataFrame(cytochrome_data)
    ref = df["Sequence"].iloc[0]
    st.dataframe(df[["Species", "Sequence"]])
    similarities = []

    for i in range(1, len(df)):
        species = df["Species"].iloc[i]
        seq = df["Sequence"].iloc[i]
        score = needleman_wunsch(ref, seq)
        a1, a2 = traceback(ref, seq, score)
        sim = sum(1 for x, y in zip(a1, a2) if x == y) / len(a1)
        similarities.append(sim)

        st.markdown(f"#### 🔎 {species} vs Human")
        st.code(f"{a1}\n{a2}\nSimilarity: {sim:.2%}")

        mutation_sites = [i for i, (x, y) in enumerate(zip(a1, a2)) if x != y]
        fig, ax = plt.subplots(figsize=(10, 1))
        ax.plot(range(len(a1)), [1] * len(a1), "|", color="gray", markersize=15)
        ax.plot(mutation_sites, [1] * len(mutation_sites), "|", color="crimson", markersize=15)
        ax.set_yticks([])
        ax.set_xlim(0, len(a1))
        ax.set_title(f"Mutation Positions: {species}")
        st.pyplot(fig)

    fig2, ax2 = plt.subplots()
    ax2.plot(df["Species"][1:], similarities, marker="o", linestyle="--", color="navy")
    ax2.set_title("Similarity Trend between Primates and Humans")
    ax2.set_ylabel("Similarity")
    ax2.set_xlabel("Species")
    ax2.grid(True)
    st.pyplot(fig2)

    st.markdown("### 🧬 Pairwise Similarity Matrix")
    matrix = []
    for seq1 in df["Sequence"]:
        row = []
        for seq2 in df["Sequence"]:
            score = needleman_wunsch(seq1, seq2)
            a1, a2 = traceback(seq1, seq2, score)
            sim = sum(x == y for x, y in zip(a1, a2)) / len(a1)
            row.append(sim)
        matrix.append(row)

    sim_df = pd.DataFrame(matrix, columns=df["Species"], index=df["Species"])
    fig3, ax3 = plt.subplots(figsize=(6, 5))
    sns.heatmap(sim_df, annot=True, cmap="viridis", fmt=".2f", square=True)
    ax3.set_title("Cytochrome C Similarity Matrix")
    st.pyplot(fig3)

    st.markdown("### 🧠 Biological Analysis & Evolutionary Inference")
    st.info("""**Why it matters:**
- Cytochrome C is extremely conserved — small changes hint at long-term evolutionary shifts.
- Closest relatives (chimpanzees, gorillas) have near-identical sequences.
- Differences occur in flexible or non-critical protein regions.
""")

# ------------------ Tab 2: Mutation ------------------
with mutation_tabs[1]:
    st.header("🧪 Mutation Case Studies")

    mutation_cases = {
        "BRCA1: c.5123C>A (p.Ala1708Glu)": {
            "ref": "AGCTGCTCAGGTCAGTACGTGACCTGAGCAGTTCAG",
            "mut": "AGCTGCTCAGGTCAGTACGAGCCTGAGCAGTTCAG",
            "info": """
**Mutation Type:** Missense (c.5123C>A → p.Ala1708Glu)  
**Effect:** Alanine → Glutamic Acid, altering folding.  
**Gene:** BRCA1 - DNA repair.  
**Impact:** High risk of breast/ovarian cancer.
"""
        },
        "BRCA1: c.68_69delAG (p.Glu23fs)": {
            "ref": "AGCTGCTCAGGTCAGTACGTGACCTGAGCAGTTCAG",
            "mut": "ACTGCTCAGGTCAGTACGTGACCTGAGCAGTTCAG",
            "info": """
**Mutation Type:** Frameshift Deletion  
**Effect:** Premature stop codon → truncated protein.  
**Gene:** BRCA1 - genome stability.  
**Impact:** Highly pathogenic.
"""
        },
        "TP53 (p.Arg248His)": {
            "ref": "MAQKTSPSGSDLGGQSLTLWKLVPRLILWQ",
            "mut": "MAQKTSPSGSDLGGQSLTLWKLVPHLILWQ",
            "info": """
**Mutation Type:** Missense  
**Effect:** Arginine → Histidine in DNA-binding region.  
**Gene:** TP53 - cell cycle regulation.  
**Impact:** Present in >50% of cancers.
"""
        },
        "COVID Spike (N501Y)": {
            "ref": "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVY",
            "mut": "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSYTRGVY",
            "info": """
**Mutation Type:** Missense  
**Effect:** N → Y in receptor-binding domain.  
**Gene:** SARS-CoV-2 spike.  
**Impact:** Increased transmission, immune escape.
"""
        }
    }

    selected = st.selectbox("🗂️ Select Case Study", list(mutation_cases.keys()))
    ref = mutation_cases[selected]["ref"]
    mut = mutation_cases[selected]["mut"]
    score = needleman_wunsch(ref, mut)
    a1, a2 = traceback(ref, mut, score)

    st.subheader("🔍 Sequence Alignment")
    st.markdown(f"<b>REF:</b> {a1}", unsafe_allow_html=True)
    st.markdown(f"<b>MUT:</b> {a2}", unsafe_allow_html=True)
    st.markdown(highlight_inline_differences(a1, a2), unsafe_allow_html=True)

    similarity = sum(r == m for r, m in zip(a1, a2)) / len(a1) * 100
    st.markdown(f"<div style='margin-top:10px'><b>Similarity:</b> {similarity:.2f}%</div>", unsafe_allow_html=True)

    st.subheader("🧬 Interpretation")
    st.info(mutation_cases[selected]["info"])
