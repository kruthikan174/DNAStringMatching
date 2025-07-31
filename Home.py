import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import time

# --- Page Config ---
st.set_page_config(page_title="DNA Aligner", layout="wide")

# --- Custom Styling ---
st.markdown("""
<style>
body {
    background-color: #f7fafc;
    font-family: 'Segoe UI', sans-serif;
    color: #333;
}

section.main > div {
    padding: 2rem;
}

h1, h2, h3 {
    color: #1e3a8a;
}

.stButton>button {
    background-color: #1d4ed8;
    color: white;
    border: none;
    border-radius: 8px;
    padding: 10px 20px;
    font-weight: 600;
    transition: 0.3s;
}

.stButton>button:hover {
    background-color: #1e40af;
    transform: scale(1.03);
}

[data-testid="stSidebar"] {
    background-color: #e0f2fe;
}

.stTextInput>div>div>input {
    border-radius: 6px;
}

hr {
    margin-top: 2rem;
    margin-bottom: 2rem;
    border: none;
    border-top: 1px solid #d1d5db;
}
</style>
""", unsafe_allow_html=True)

# --- Sidebar: Scoring Settings ---
st.sidebar.header("⚙️ Settings")
mode = st.sidebar.radio("Input Mode", ["Manual", "Upload CSV"], help="Choose how to input DNA sequences")
match = st.sidebar.number_input("Match Score", value=1, help="Score for matching bases")
mismatch = st.sidebar.number_input("Mismatch Penalty", value=-1, help="Penalty for mismatched bases")
gap = st.sidebar.number_input("Gap Penalty", value=-2, help="Penalty for gaps")

# --- Title & Intro ---
st.title("🧬 DNA Sequence Alignment Tool")
st.markdown("""
Welcome to the **interactive sequence alignment dashboard**. This tool helps you:
- Compare sequences globally and locally
- Visualize scoring matrices
- Highlight mismatches and similarity

Choose input, adjust scoring, and click **Run Alignment**!
""")

# --- Algorithm Overview ---
with st.expander("📚 About Algorithms", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### 🔵 Needleman-Wunsch")
        st.info("""
- **Global alignment** of full sequences
- Suitable for complete similarity comparison
""")
    with col2:
        st.markdown("### 🟢 Smith-Waterman")
        st.success("""
- **Local alignment** for best matching
- Ideal for short matches
""")

# --- DNA Validator ---
def is_valid_dna(seq):
    return bool(re.fullmatch(r"[ATCG]+", seq.upper()))

# --- Alignment Algorithms ---
def needleman_wunsch(s1, s2, match, mismatch, gap):
    n, m = len(s1), len(s2)
    dp = np.zeros((n+1, m+1), dtype=int)
    for i in range(n+1): dp[i][0] = i * gap
    for j in range(m+1): dp[0][j] = j * gap
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = dp[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch)
            dp[i][j] = max(diag, dp[i-1][j] + gap, dp[i][j-1] + gap)
    a1, a2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        score = dp[i][j]
        if i > 0 and j > 0 and score == dp[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch):
            a1, a2 = s1[i-1] + a1, s2[j-1] + a2; i -= 1; j -= 1
        elif i > 0 and score == dp[i-1][j] + gap:
            a1, a2 = s1[i-1] + a1, "-" + a2; i -= 1
        else:
            a1, a2 = "-" + a1, s2[j-1] + a2; j -= 1
    return dp, a1, a2, dp[n][m]

def smith_waterman(s1, s2, match, mismatch, gap):
    n, m = len(s1), len(s2)
    dp = np.zeros((n+1, m+1), dtype=int)
    max_score, max_pos = 0, (0, 0)
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = dp[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch)
            dp[i][j] = max(0, diag, dp[i-1][j] + gap, dp[i][j-1] + gap)
            if dp[i][j] > max_score:
                max_score, max_pos = dp[i][j], (i, j)
    i, j = max_pos
    a1, a2 = "", ""
    while i > 0 and j > 0 and dp[i][j] > 0:
        score = dp[i][j]
        if score == dp[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch):
            a1, a2 = s1[i-1] + a1, s2[j-1] + a2; i -= 1; j -= 1
        elif score == dp[i-1][j] + gap:
            a1, a2 = s1[i-1] + a1, "-" + a2; i -= 1
        else:
            a1, a2 = "-" + a1, s2[j-1] + a2; j -= 1
    return dp, a1, a2, max_score

# --- Mismatch Highlighter ---
def highlight_mismatches(a1, a2):
    out1, out2 = "", ""
    for x, y in zip(a1, a2):
        if x != y:
            out1 += f":red[{x}]"
            out2 += f":red[{y}]"
        else:
            out1 += x
            out2 += y
    return out1, out2

# --- Input Section ---
st.markdown("### 📥 DNA Input")

seq_pairs = []
if mode == "Manual":
    s1 = st.text_input("🔡 Sequence 1").upper()
    s2 = st.text_input("🔡 Sequence 2").upper()
    if s1 and s2:
        if is_valid_dna(s1) and is_valid_dna(s2):
            seq_pairs = [(s1, s2)]
        else:
            st.error("❌ Invalid characters! Use A, T, C, G only.")
else:
    file = st.file_uploader("📂 Upload CSV (2 columns: Sequence1, Sequence2)", type="csv")
    if file:
        df = pd.read_csv(file, header=None)
        for row in df.itertuples(index=False):
            if is_valid_dna(row[0]) and is_valid_dna(row[1]):
                seq_pairs.append((row[0], row[1]))
        st.success(f"✅ Loaded {len(seq_pairs)} valid sequence pairs.")

# --- Alignment ---
if st.button("🚀 Run Alignment") and seq_pairs:
    for idx, (s1, s2) in enumerate(seq_pairs, start=1):
        st.markdown(f"## 🧪 Results for Pair {idx}")
        col1, col2 = st.columns(2)

        with st.spinner("Aligning sequences..."):
            time.sleep(0.6)

        # Needleman-Wunsch
        with col1:
            st.markdown("#### 🔵 Global Alignment")
            t0 = time.time()
            m, a1, a2, score = needleman_wunsch(s1, s2, match, mismatch, gap)
            dur = time.time() - t0
            matches = sum(1 for x, y in zip(a1, a2) if x == y)
            sim = matches / len(s1)
            h1, h2 = highlight_mismatches(a1, a2)
            st.markdown(f"**Seq1:** {h1}")
            st.markdown(f"**Seq2:** {h2}")
            st.info(f"Score: `{score}`, Matches: `{matches}`, Similarity: `{sim:.2%}`")
            with st.expander("📊 Score Matrix"):
                fig, ax = plt.subplots(figsize=(5, 4))
                sns.heatmap(m, annot=True, fmt="d", cmap="Blues",
                            xticklabels=["-"] + list(s2), yticklabels=["-"] + list(s1), ax=ax)
                st.pyplot(fig)

        # Smith-Waterman
        with col2:
            st.markdown("#### 🟢 Local Alignment")
            t1 = time.time()
            m, a1, a2, score = smith_waterman(s1, s2, match, mismatch, gap)
            dur = time.time() - t1
            matches = sum(1 for x, y in zip(a1, a2) if x == y)
            sim = matches / len(s1)
            h1, h2 = highlight_mismatches(a1, a2)
            st.markdown(f"**Seq1:** {h1}")
            st.markdown(f"**Seq2:** {h2}")
            st.success(f"Score: `{score}`, Matches: `{matches}`, Similarity: `{sim:.2%}`")
            with st.expander("📊 Score Matrix"):
                fig, ax = plt.subplots(figsize=(5, 4))
                sns.heatmap(m, annot=True, fmt="d", cmap="Greens",
                            xticklabels=["-"] + list(s2), yticklabels=["-"] + list(s1), ax=ax)
                st.pyplot(fig)
else:
    st.info("📌 Please provide sequences and click **Run Alignment**.")
