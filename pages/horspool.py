import streamlit as st
import re

# Scoring
def local_alignment_horspool_like(text, pattern, match_score=1, mismatch_penalty=-1):
    text = text.upper()
    pattern = pattern.upper()
    best_score = float("-inf")
    results = []

    for i in range(len(text) - len(pattern) + 1):
        sub = text[i:i+len(pattern)]
        score = 0
        for a, b in zip(sub, pattern):
            if a == b:
                score += match_score
            else:
                score += mismatch_penalty
        results.append((i, sub, pattern, score))
        best_score = max(best_score, score)

    results.sort(key=lambda x: x[3], reverse=True)
    return results[:5], best_score

# Coloring Matches
def highlight(a1, a2):
    s1, s2, match = "", "", ""
    for x, y in zip(a1, a2):
        if x == y:
            s1 += f":green[{x}]"
            s2 += f":green[{y}]"
            match += "|"
        else:
            s1 += f":red[{x}]"
            s2 += f":red[{y}]"
            match += "."
    return s1, s2, match

# Streamlit UI
st.set_page_config("Horspool-like Local Alignment", layout="centered")
st.title("Horspool: Local DNA Alignment")
text = st.text_input("DNA Text (Long Sequence)").upper()
pattern = st.text_input("Pattern to Match (Short Sequence)").upper()

match_score = 2
mismatch_penalty = -1

if st.button("Run Alignment"):
    if not re.fullmatch(r"[ATCG]+", text) or not re.fullmatch(r"[ATCG]+", pattern):
        st.error("Please use only A, T, C, G in sequences.")
    elif len(pattern) > len(text):
        st.warning("Pattern is longer than the text.")
    else:
        top_hits, best = local_alignment_horspool_like(text, pattern, match_score, mismatch_penalty)
        st.markdown("### Top Matches")
        for i, (pos, sub, pat, score) in enumerate(top_hits, 1):
            h1, h2, match_line = highlight(sub, pat)
            st.markdown(f"**{i}. Position `{pos}` | Score: `{score}`**")
            st.markdown(f"**Text  :** {h1}")
            st.markdown(f"**Pattern:** {h2}")
            st.markdown("---")
