🧬 Bioinformatics Sequence Analysis Toolkit

This repository contains an interactive Streamlit-based toolkit for exploring biological sequence analysis.
It combines classical alignment algorithms, mutation simulation, and case studies to help students and researchers understand molecular evolution.

📌 Features

🔹 1. Sequence Alignment (Home.py)

Global Alignment (Needleman–Wunsch) – Aligns entire sequences with gap penalties.

Local Alignment (Smith–Waterman) – Finds the best-matching subsequences.

Customizable Parameters – Match, mismatch, and gap penalties.

Outputs:

Alignment score

Aligned sequences with mismatches highlighted

Dynamic programming (DP) score matrix heatmap

🔹 2. Horspool Pattern Search (horspool.py)

Efficient pattern searching inspired by Horspool’s algorithm.

Compares substrings between two sequences for local alignment.

Highlights mismatches, gap penalties, and best matching regions.

Visualizations:

Alignment score heatmap

Score distribution histogram

🔹 3. Mutation Simulation Tools (Mutation Simulation Tools.py)

Manual mutation editor – Introduce substitutions, insertions, deletions.

Random mutation generator – Simulates mutation probability across sequences.

Alignment comparison – See how mutations affect global alignment score.

Visualizations:

Mutation effects on alignment

Distribution of random mutation scores

🔹 4. Case Studies: Cytochrome C Evolution

Compare Cytochrome C protein sequences across primates (human, chimp, gorilla, orangutan, macaque).

Analyze evolutionary relationships using alignment scores.

Visual outputs:

Heatmap of pairwise alignment scores

Bar plot of sequence similarity

📊 Tech Stack

Python – Core implementation

Streamlit – Interactive web-based UI

NumPy / Pandas – Data handling & calculations

Matplotlib / Seaborn – Heatmaps & plots

Biopython (optional future extension) – Advanced sequence utilities

🚀 Future Enhancements

Add frameshift & codon-level mutations.

Support for protein translation (DNA → amino acids).

Expand case studies (Hemoglobin mutations, viral genome comparisons).

Optimize alignment algorithms with NumPy/Numba for large sequences.

📖 Educational Use

This toolkit is designed for:

Bioinformatics students learning alignment algorithms.

Researchers experimenting with sequence analysis.

Educators demonstrating molecular evolution concepts.
