import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import re

# ✅ Set page config must be first Streamlit command
st.set_page_config(page_title="BioAlign - Sequence Alignment Tool", layout="wide")

# ----------------- Top-Level Tabs -----------------
main_tabs = st.tabs(["🧬 DNA Matching", "🧬 Mutation"])

# ----------------- DNA Matching -----------------
with main_tabs[0]:
    # ✅ Background styling
    st.markdown(
        """
        <style>
        div[data-testid="stApp"] {
            background-image: url("https://t3.ftcdn.net/jpg/04/35/46/22/360_F_435462259_0t4ktQ8K1oNSxbUbkJpYrkIO9xOFB7t2.jpg");
            background-size: cover;
            background-repeat: no-repeat;
            background-attachment: fixed;
            background-position: center;
            color: white;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

    st.title("🔬 BioALign - DNA Sequence Alignment Tool")
    st.markdown("Compare sequences using **Needleman-Wunsch** or **Smith-Waterman** algorithms.")

    # Alignment algorithms
    def needleman_wunsch(seq1, seq2, match, mismatch, gap):
        n, m = len(seq1), len(seq2)
        score = np.zeros((n + 1, m + 1), dtype=int)
        for i in range(n + 1): score[i][0] = gap * i
        for j in range(m + 1): score[0][j] = gap * j
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
                delete = score[i - 1][j] + gap
                insert = score[i][j - 1] + gap
                score[i][j] = max(diag, delete, insert)
        align1, align2 = "", ""
        i, j = n, m
        while i > 0 or j > 0:
            current = score[i][j]
            if i > 0 and j > 0 and current == score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
                align1 = seq1[i - 1] + align1
                align2 = seq2[j - 1] + align2
                i -= 1
                j -= 1
            elif i > 0 and current == score[i - 1][j] + gap:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1
        return score, align1, align2, score[n][m]

    def smith_waterman(seq1, seq2, match, mismatch, gap):
        n, m = len(seq1), len(seq2)
        score = np.zeros((n + 1, m + 1), dtype=int)
        max_score, max_pos = 0, (0, 0)
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
                delete = score[i - 1][j] + gap
                insert = score[i][j - 1] + gap
                score[i][j] = max(0, diag, delete, insert)
                if score[i][j] > max_score:
                    max_score = score[i][j]
                    max_pos = (i, j)
        align1, align2 = "", ""
        i, j = max_pos
        while i > 0 and j > 0 and score[i][j] != 0:
            current = score[i][j]
            if current == score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
                align1 = seq1[i - 1] + align1
                align2 = seq2[j - 1] + align2
                i -= 1
                j -= 1
            elif current == score[i - 1][j] + gap:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1
        return score, align1, align2, max_score

    def is_valid_dna(seq):
        return bool(re.fullmatch(r'[ATCG]+', seq.upper()))

    def read_multiple_sequences(file):
        df = pd.read_csv(file, header=None)
        sequences = df.apply(lambda row: ''.join(row.dropna().astype(str)).strip().upper(), axis=1)
        return [s for s in sequences if is_valid_dna(s)]

    with st.sidebar:
        st.header("🔧 Parameters")
        algorithm = st.radio("Select Algorithm", ["Needleman-Wunsch", "Smith-Waterman"])
        match = st.number_input("Match Score", value=1)
        mismatch = st.number_input("Mismatch Penalty", value=-1)
        gap = st.number_input("Gap Penalty", value=-2)
        upload_option = st.radio("Input Mode", ["Manual Entry", "Upload CSV Files"])

    seq_list1, seq_list2 = [], []
    if upload_option == "Manual Entry":
        seq1 = st.text_input("Enter Sequence 1 (e.g. AGCT)", "").upper()
        seq2 = st.text_input("Enter Sequence 2 (e.g. AGTT)", "").upper()
        if seq1 and seq2:
            if not is_valid_dna(seq1) or not is_valid_dna(seq2):
                st.warning("Only characters A, T, C, G are allowed in DNA sequences.")
            else:
                seq_list1 = [seq1]
                seq_list2 = [seq2]
    else:
        file1 = st.file_uploader("Upload CSV for Sequence Set 1", type="csv")
        file2 = st.file_uploader("Upload CSV for Sequence Set 2", type="csv")
        if file1 and file2:
            try:
                seq_list1 = read_multiple_sequences(file1)
                seq_list2 = read_multiple_sequences(file2)
                if len(seq_list1) != len(seq_list2):
                    st.error("CSV files must contain the same number of valid sequences.")
                    seq_list1, seq_list2 = [], []
                elif len(seq_list1) == 0:
                    st.error("No valid DNA sequences found in the uploaded files.")
                else:
                    st.success(f"✅ {len(seq_list1)} valid sequence pairs loaded!")
            except Exception as e:
                st.error(f"Error reading files: {e}")

    if st.button("🚀 Run Alignment"):
        if seq_list1 and seq_list2:
            for idx, (s1, s2) in enumerate(zip(seq_list1, seq_list2), start=1):
                st.markdown(f"### 🔬 Pair {idx}")
                start_time = time.time()
                matrix, a1, a2, score = (
                    needleman_wunsch(s1, s2, match, mismatch, gap)
                    if algorithm.startswith("Needleman")
                    else smith_waterman(s1, s2, match, mismatch, gap)
                )
                exec_time = time.time() - start_time
                matches = sum(1 for x, y in zip(a1, a2) if x == y)
                alignment_length = len(a1)
                st.code(f"Seq1: {a1}\nSeq2: {a2}\nScore: {score}\nMatches: {matches}\nLength: {alignment_length}\nTime: {exec_time:.4f}s")
                with st.expander("📊 Show Matrix"):
                    fig, ax = plt.subplots(figsize=(6, 4))
                    sns.heatmap(matrix, annot=False, fmt="d", cmap="coolwarm", ax=ax,
                                xticklabels=["-"] + list(s2),
                                yticklabels=["-"] + list(s1))
                    ax.set_title(f"Score Matrix for Pair {idx}")
                    st.pyplot(fig)
                similarity = matches / alignment_length if alignment_length > 0 else 0
                if similarity > 0.7:
                    st.success("🟢 Highly Similar")
                elif similarity > 0.4:
                    st.warning("🟠 Moderately Similar")
                else:
                    st.error("🔴 Low Similarity")
        else:
            st.warning("Please enter or upload valid sequences to proceed.")

# ----------------- Mutation Simulators -----------------
with main_tabs[1]:
    st.subheader("🧬 Mutation Simulation Tools")
    mutation_tabs = st.tabs(["🔬 Error Injection", "🔄 Sequence Evolution", "📚 Case Studies"])

    # Error Injection
    with mutation_tabs[0]:
        st.header("🧪 Error Injection Simulator")
        original_seq = st.text_input("Original DNA or Protein Sequence").upper()
        mutation_type = st.selectbox("Mutation Type", ["Insertion", "Deletion", "Substitution"])
        position = st.number_input("Position (0-based)", min_value=0, value=0)
        base = st.text_input("New Base (for Insertion/Substitution)").upper()
        inject = st.button("Apply Mutation")

        def inject_mut(seq, mut_type, pos, base):
            if pos < 0 or pos > len(seq): return "Invalid position", ""
            if mut_type == "Insertion": return "", seq[:pos] + base + seq[pos:]
            if mut_type == "Deletion": return "", seq[:pos] + seq[pos + 1:]
            if mut_type == "Substitution": return "", seq[:pos] + base + seq[pos + 1:]
            return "Unknown mutation", ""

        if inject:
            err, mutated = inject_mut(original_seq, mutation_type, position, base)
            if err:
                st.error(err)
            else:
                st.code(f"Original: {original_seq}\nMutated : {mutated}")
                _, a1, a2, score = needleman_wunsch(original_seq, mutated, match, mismatch, gap)
                st.code(f"Alignment:\n{a1}\n{a2}\nScore: {score}")

    # Sequence Evolution
    with mutation_tabs[1]:
        st.header("🔄 Sequence Evolution Explorer")

        sub_tabs = st.tabs(["🧪 Simulate Mutation", "🌍Evolution Case Study"])

        # --- Existing Simulation Code Here ---
        with sub_tabs[0]:
            st.header("🔄 Sequence Evolution Simulator")
            base_seq = st.text_input("Base Sequence").upper()
            mutation_pct = st.slider("Mutation %", 0, 100, 10)
            simulate = st.button("Simulate Evolution")


            def simulate_evolution(seq, percent):
                seq = list(seq)
                total = len(seq)
                changes = int((percent / 100) * total)
                for _ in range(changes):
                    i = np.random.randint(0, total)
                    seq[i] = np.random.choice([b for b in "ATCG" if b != seq[i]])
                return ''.join(seq)


            if simulate:
                if not base_seq:
                    st.warning("Enter a valid base sequence")
                else:
                    evolved = simulate_evolution(base_seq, mutation_pct)
                    st.subheader("🔁 Evolution Summary")
                    st.code(f"Original: {base_seq}\nEvolved : {evolved}")

                    # Alignment
                    _, a1, a2, score = needleman_wunsch(base_seq, evolved, match, mismatch, gap)
                    matches = sum(1 for x, y in zip(a1, a2) if x == y)
                    similarity = matches / len(a1)

                    st.subheader("📊 Alignment & Mutation Impact")
                    st.code(f"{a1}\n{a2}")
                    st.metric("Alignment Score", score)
                    st.metric("Similarity", f"{similarity:.2f}")

                    # Feedback
                    if similarity > 0.7:
                        st.success("🟢 Minimal impact. Sequence still closely resembles the original.")
                    elif similarity > 0.4:
                        st.warning("🟠 Moderate divergence. Potential functional impact.")
                    else:
                        st.error("🔴 High divergence. Functionality may be lost or drastically altered.")

                    # Match/mismatch chart
                    fig, ax = plt.subplots()
                    ax.bar(["Match", "Mismatch"], [matches, len(a1) - matches], color=["green", "crimson"])
                    ax.set_title("Base-by-Base Comparison")
                    st.pyplot(fig)

                    # Mutation % vs similarity plot
                    st.subheader("📈 Mutation % vs Similarity")
                    x_vals, y_vals = [], []
                    for p in range(0, 101, 10):
                        m_seq = simulate_evolution(base_seq, p)
                        _, a1, a2, _ = needleman_wunsch(base_seq, m_seq, match, mismatch, gap)
                        sim = sum(x == y for x, y in zip(a1, a2)) / len(a1)
                        x_vals.append(p)
                        y_vals.append(sim)
                    fig2, ax2 = plt.subplots()
                    ax2.plot(x_vals, y_vals, marker='o', linestyle='-', color='blue')
                    ax2.set_title("Similarity vs Mutation %")
                    ax2.set_xlabel("Mutation %")
                    ax2.set_ylabel("Similarity")
                    ax2.grid(True)
                    st.pyplot(fig2)

                    # Explanation
                    st.info("""
                    🔍 **Interpretation:** As the mutation percentage increases, the similarity typically decreases.  
                    This helps visualize how random mutations can erode genetic similarity and affect biological function.
                    """)
        # --- Real Evolution Study ---
        with sub_tabs[1]:
            st.subheader("🌍 Evolution of Primates via Cytochrome C Protein")

            # Static dataset
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
            ref = df["Sequence"].iloc[0]  # Human as reference

            st.dataframe(df[["Species", "Sequence"]])

            similarities = []
            mutation_maps = []

            st.markdown("### 🧬 Stepwise Comparison to Human Sequence")

            for i in range(1, len(df)):
                species = df["Species"].iloc[i]
                seq = df["Sequence"].iloc[i]
                _, a1, a2, score = needleman_wunsch(ref, seq, match, mismatch, gap)
                sim = sum(1 for x, y in zip(a1, a2) if x == y) / len(a1)
                similarities.append(sim)

                st.markdown(f"#### 🔎 {species} vs Human")
                st.code(f"{a1}\n{a2}\nSimilarity: {sim:.2f} | Alignment Score: {score}")

                mutation_sites = [i for i, (x, y) in enumerate(zip(a1, a2)) if x != y]
                mutation_maps.append(mutation_sites)

                # Mutation position map
                fig, ax = plt.subplots(figsize=(10, 1))
                ax.plot(range(len(a1)), [1] * len(a1), "|", color="gray", markersize=15, label="Reference")
                ax.plot(mutation_sites, [1] * len(mutation_sites), "|", color="crimson", markersize=15,
                        label="Mutations")
                ax.set_yticks([])
                ax.set_xlim(0, len(a1))
                ax.set_title(f"Mutation Positions: {species}")
                ax.legend()
                st.pyplot(fig)

            # Evolution similarity plot
            st.markdown("### 📈 Evolutionary Similarity to Human")
            fig2, ax2 = plt.subplots()
            ax2.plot(df["Species"][1:], similarities, marker="o", linestyle="--", color="navy")
            ax2.set_title("Similarity Trend from Primates to Human")
            ax2.set_ylabel("Similarity")
            ax2.set_xlabel("Species")
            ax2.grid(True)
            st.pyplot(fig2)

            # Heatmap of all pairwise similarities
            st.markdown("### 🧬 Pairwise Similarity Matrix")
            matrix = []
            for seq1 in df["Sequence"]:
                row = []
                for seq2 in df["Sequence"]:
                    _, a1, a2, _ = needleman_wunsch(seq1, seq2, match, mismatch, gap)
                    sim = sum(x == y for x, y in zip(a1, a2)) / len(a1)
                    row.append(sim)
                matrix.append(row)

            sim_df = pd.DataFrame(matrix, columns=df["Species"], index=df["Species"])
            fig3, ax3 = plt.subplots(figsize=(6, 5))
            sns.heatmap(sim_df, annot=True, cmap="viridis", fmt=".2f", square=True)
            ax3.set_title("Cytochrome C Similarity Matrix")
            st.pyplot(fig3)
        # --- Interpretation Section ---
            st.markdown("### 🧠 Biological Analysis & Evolutionary Inference")
            st.info("""
                    ### 🔬 What Do These Results Mean?

                    - **Cytochrome C** is a highly conserved protein involved in cellular respiration.
                    - Despite millions of years of divergence, primates show **>95% similarity** in Cytochrome C, indicating that **mutations in this gene are rare and often deleterious**.

                    ---

                    ### 📌 Observations:

                    - **Chimpanzees and Gorillas** have the highest similarity to humans (>98%), suggesting **very recent evolutionary divergence**.
                    - **Orangutans and Macaques** show slightly more mutations — these correlate with their earlier divergence in the primate tree.
                    - The **mutation position maps** indicate that most substitutions occur in **non-active** regions of the protein, minimizing functional disruption.

                    ---

                    ### 🧬 Mutation Analysis:

                    - Most detected mutations are **conservative** (amino acids with similar properties), minimizing functional impact.
                    - **Rare non-conservative changes** (e.g., polar ➝ nonpolar) could reflect **adaptations or neutral drift** over time.
                    - The fact that **core regions are identical across species** reinforces their critical role in protein function.

                    ---

                    ### 🧭 Evolutionary Insight:

                    This analysis reinforces the concept of **molecular clock** — the idea that genes mutate at relatively steady rates over time.  
                    Highly conserved genes like Cytochrome C allow us to:
                    - Estimate **evolutionary distances**
                    - Identify **functionally critical regions**
                    - Understand **selective pressure** across species

                    In short: **Few changes = essential gene**.  
                    More mutations? It may be **less constrained** or **species-adaptive**.

                    """)


    # Real Dataset Case Studies
    with mutation_tabs[2]:
        st.header("📚 Case Studies")
        choice = st.selectbox("Dataset Type", ["COVID Spike Protein", "TP53 Protein Mutations"])
        file = st.file_uploader("Upload CSV", type="csv", key=choice)

        if file:
            df = pd.read_csv(file)
            st.dataframe(df)

            similarities = []
            labels = []

            for index, row in df.iterrows():
                st.subheader(f"🧬 {row.get('Mutation', row.get('Variant', f'Case {index + 1}'))}")
                st.caption(row.get('Description', ''))

                ref = row['Reference'].upper()
                mut = row['Mutated'].upper()
                implication = row.get("Implication", "")
                risk_level = row.get("Risk_Level", "Unknown")

                # Alignment and similarity
                _, a1, a2, score = needleman_wunsch(ref, mut, match, mismatch, gap)
                matches = sum(1 for x, y in zip(a1, a2) if x == y)
                similarity = matches / len(a1)

                st.code(f"{a1}\n{a2}")
                st.metric("Alignment Score", score)
                st.metric("Similarity", f"{similarity:.2f}")

                # Scientific interpretation
                st.markdown("### 🧠 Biological & Clinical Impact")

                if implication:
                    st.info(f"**Implication:** {implication}")
                else:
                    if "TP53" in choice and "R175H" in row.get("Mutation", ""):
                        st.warning(
                            "🧬 This hotspot TP53 mutation is associated with structural loss in the DNA-binding domain and linked to several cancers.")
                    elif "D614G" in row.get("Mutation", ""):
                        st.info(
                            "🦠 The D614G mutation in the SARS-CoV-2 spike protein is associated with enhanced transmissibility and infectivity.")
                    else:
                        st.write("📌 No direct clinical annotation available for this mutation.")

                st.markdown(f"**🧪 Functional Classification:** `{risk_level}`")
                if risk_level.lower() == "high":
                    st.error("⚠️ Likely to significantly disrupt protein function or contribute to disease.")
                elif risk_level.lower() == "moderate":
                    st.warning("🟠 May impair protein function or structure moderately.")
                else:
                    st.success("🟢 Mutation is likely benign or low impact.")

                # Mutation position map
                st.markdown("**📍 Mutation Position Map**")
                mutation_sites = [i for i, (r, m) in enumerate(zip(ref, mut)) if r != m]
                fig, ax = plt.subplots(figsize=(10, 1))
                ax.plot(range(len(ref)), [1] * len(ref), "|", color="lightgray", markersize=15, label="Reference")
                ax.plot(mutation_sites, [1] * len(mutation_sites), "|", color="crimson", markersize=15, label="Mutated")
                ax.set_yticks([])
                ax.set_xlim(0, len(ref))
                ax.set_title("Mutation Positions")
                ax.legend()
                st.pyplot(fig)

                # Amino acid property shift (only if protein)
                if all(c in "ACDEFGHIKLMNPQRSTVWY" for c in ref + mut):
                    aa_props = {
                        'A': 'nonpolar', 'C': 'polar', 'D': 'acidic', 'E': 'acidic', 'F': 'nonpolar',
                        'G': 'nonpolar', 'H': 'basic', 'I': 'nonpolar', 'K': 'basic', 'L': 'nonpolar',
                        'M': 'nonpolar', 'N': 'polar', 'P': 'nonpolar', 'Q': 'polar', 'R': 'basic',
                        'S': 'polar', 'T': 'polar', 'V': 'nonpolar', 'W': 'nonpolar', 'Y': 'polar'
                    }
                    ref_props = [aa_props.get(c, 'unknown') for c in ref]
                    mut_props = [aa_props.get(c, 'unknown') for c in mut]
                    change = [int(r != m) for r, m in zip(ref_props, mut_props)]
                    fig2, ax2 = plt.subplots(figsize=(10, 2))
                    ax2.plot(change, color='purple')
                    ax2.set_title("Amino Acid Property Change Map")
                    ax2.set_xlabel("Position")
                    ax2.set_yticks([0, 1])
                    ax2.set_yticklabels(["Same", "Changed"])
                    st.pyplot(fig2)

                # Collect data for heatmap
                similarities.append(similarity)
                labels.append(row.get('Mutation', f"Case {index + 1}"))

            # Similarity matrix
            st.markdown("### 🧬 Overall Similarity Heatmap")
            sim_df = pd.DataFrame([similarities], columns=labels)
            fig3, ax3 = plt.subplots(figsize=(len(labels) * 0.6, 1.5))
            sns.heatmap(sim_df, cmap="YlGnBu", annot=True, fmt=".2f", cbar=False,
                        xticklabels=labels, yticklabels=["Similarity"])
            ax3.set_title("Mutation-to-Reference Similarity")
            st.pyplot(fig3)

