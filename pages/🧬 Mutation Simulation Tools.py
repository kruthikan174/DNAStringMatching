import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =================== Needleman-Wunsch Alignment ===================
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    score = [[0] * (n + 1) for _ in range(m + 1)]
    pointer = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        score[i][0] = gap * i
    for j in range(n + 1):
        score[0][j] = gap * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = score[i - 1][j] + gap
            insert = score[i][j - 1] + gap
            score[i][j] = max(diag, delete, insert)
            pointer[i][j] = [diag, delete, insert].index(score[i][j]) + 1

    aligned1, aligned2 = '', ''
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and pointer[i][j] == 1:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif i > 0 and(j == 0 or pointer[i][j] == 2):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = '-' + aligned2
            i -= 1
        elif j > 0:
            aligned1 = '-' + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1

    return score[m][n], aligned1, aligned2

# =================== Streamlit Page Setup ===================
match, mismatch, gap = 1, -1, -2
st.set_page_config(page_title="Mutation Tools", layout="wide")

st.title("🧬 Mutation Simulation Tools")
st.markdown("""
Welcome to the **Mutation Simulation Lab**, where you can:

- 🔄 Inject genetic errors (Insertion, Deletion, Substitution)
- 🔬 Observe mutation impact via alignment
- 🧪 Simulate evolution by applying random mutations
""")

# =================== Tabs Setup ===================
tabs = st.tabs(["🔬 Error Injection", "🔄 Sequence Evolution"])

# =================== Error Injection Tab ===================
with tabs[0]:
    st.header("🧪 Manual Mutation Injection")
    st.markdown("""
Simulate specific mutations by selecting the type of error and observing how the sequence changes.

- **Insertion** adds a new base at a given position  
- **Deletion** removes a base  
- **Substitution** replaces an existing base  
""")

    original_seq = st.text_input("🔤 Original DNA/Protein Sequence").upper()
    mutation_type = st.selectbox("⚙️ Mutation Type", ["Insertion", "Deletion", "Substitution"])

    position = st.number_input("📍 Position (0-based)", min_value=0, value=0)

    # Show base input only for Insertion and Substitution
    if mutation_type in ["Insertion", "Substitution"]:
        base = st.text_input("🧬 New Base (for Insertion/Substitution)").upper()
    else:
        base = ""  # Hide and clear base if not required

    inject = st.button("💥 Apply Mutation")

    def inject_mut(seq, mut_type, pos, base):
        if not seq:
            return "Sequence is empty", ""
        if pos < 0 or pos > len(seq):
            return "Invalid position", ""
        if mut_type == "Insertion":
            if not base: return "Base required for insertion", ""
            return "", seq[:pos] + base + seq[pos:]
        if mut_type == "Deletion":
            if pos >= len(seq): return "Position out of range for deletion", ""
            return "", seq[:pos] + seq[pos + 1:]
        if mut_type == "Substitution":
            if pos >= len(seq): return "Position out of range for substitution", ""
            if not base: return "Base required for substitution", ""
            return "", seq[:pos] + base + seq[pos + 1:]
        return "Unknown mutation", ""

    if inject:
        err, mutated = inject_mut(original_seq, mutation_type, position, base)
        if err:
            st.error(err)
        else:
            st.code(f"Original: {original_seq}\nMutated : {mutated}")

            _, a1, a2 = needleman_wunsch(original_seq, mutated, match, mismatch, gap)
            st.subheader("🔍 Alignment After Mutation")
            st.code(f"{a1}\n{a2}")

            similarity = sum(x == y for x, y in zip(a1, a2)) / len(a1)
            st.metric("🔗 Similarity", f"{similarity:.2f}")

            st.info("""
**Why Alignment Matters:**

Even a small mutation can cause large downstream shifts in alignment, affecting protein folding or function.

- Insertions and deletions are more disruptive than substitutions.
- Observe whether mutations break long matches or introduce shifts.
""")


# =================== Evolution Tab ===================
with tabs[1]:
    st.header("🔄 Sequence Evolution Simulator")
    st.markdown("""
This tool simulates random mutations over a sequence to mimic **natural evolution**.

You can:
- Apply a specific % of random mutations
- Visualize alignment changes
- Plot how similarity drops with increasing mutation levels
""")

    base_seq = st.text_input("🌱 Base Sequence").upper()
    mutation_pct = st.slider("🎚️ Mutation Percentage", 0, 100, 10)
    simulate = st.button("🚀 Simulate Evolution")

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
            st.warning("⚠️ Please enter a valid sequence.")
        else:
            evolved = simulate_evolution(base_seq, mutation_pct)
            st.subheader("🔁 Evolution Summary")
            st.code(f"Original: {base_seq}\nEvolved : {evolved}")

            _, a1, a2 = needleman_wunsch(base_seq, evolved, match, mismatch, gap)
            matches = sum(1 for x, y in zip(a1, a2) if x == y)
            similarity = matches / len(a1)

            st.subheader("📊 Alignment & Mutation Impact")
            st.code(f"{a1}\n{a2}")
            st.metric("🔗 Similarity", f"{similarity:.2f}")
            st.metric("🧬 Matching Bases", f"{matches} / {len(a1)}")

            if similarity > 0.7:
                st.success("🟢 Minimal impact — sequence is still mostly intact.")
            elif similarity > 0.4:
                st.warning("🟠 Moderate divergence — could cause partial functional loss.")
            else:
                st.error("🔴 High divergence — likely to disrupt function.")

            fig, ax = plt.subplots()
            ax.bar(["Match", "Mismatch"], [matches, len(a1) - matches], color=["green", "crimson"])
            ax.set_title("Base-by-Base Comparison")
            st.pyplot(fig)

            st.subheader("📈 Mutation % vs Similarity Curve")
            x_vals, y_vals = [], []
            for p in range(0, 101, 5):
                m_seq = simulate_evolution(base_seq, p)
                _, a1, a2 = needleman_wunsch(base_seq, m_seq, match, mismatch, gap)
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

            st.info("""
**Interpretation:**

- The curve demonstrates how increasing mutations leads to **genetic drift**.
- Sharp drops indicate sensitive regions where small mutations break alignment.
- This reflects how genetic diversity builds up during evolution over generations.
""")
