"""Exploratory overlap of genomic-layer genes (GWAS Catalog + ClinVar) across diseases.

Harmonises gene symbols to current HGNC nomenclature, counts gene entities per
disease, tabulates genes shared by >=2 diseases (all entities and protein-coding
genes), draws Venn diagrams and runs a size-matched permutation test of the
number of genes shared by >=k diseases (no neuro/ophthalmic criterion; the
criterion-based analysis used in the manuscript is in 06_genomics_and_supp.py).

Usage (from the repository root):
    python scripts/05_genomics_overlap.py
"""
import itertools
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from venn import venn  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "data" / "genomics"
HGNC = ROOT / "data" / "reference" / "hgnc_protein_coding.txt"
OUT = ROOT / "results" / "genomics" / "exploratory"
OUT.mkdir(parents=True, exist_ok=True)

# (GWAS Catalog file, ClinVar file); DR has no ClinVar records
FILES = {
    "AD": ("ADgenesGWAS.csv", "ADgenesClinvar.csv"),
    "PD": ("PDgenesGWAS.csv", "PDGenesClinvar.csv"),
    "AMD": ("AMDgenesGWAS.csv", "ARMDgenesClinvar.csv"),
    "DR": ("DRgenesGWAS.csv", None),
    "Glaucoma": ("GlaucomagenesGWAS.csv", "GlaucomaGenesClinVar.csv"),
}
LAB = {"AD": "Alzheimer", "PD": "Parkinson", "AMD": "AMD", "DR": "Diabetic retinopathy", "Glaucoma": "Glaucoma"}

# ---- HGNC protein-coding genes and previous-symbol map ----
hg = pd.read_csv(HGNC, sep="\t", dtype=str)
official = set(hg.symbol)
prev = {}
for sym, p in zip(hg.symbol, hg.prev_symbol.fillna("")):
    for q in p.split("|"):
        if q:
            prev[q] = sym


def harmonise(x):
    """Map a previous symbol to its current HGNC symbol."""
    return prev[x] if (x not in official and x in prev) else x


def load(fname):
    if fname is None:
        return set()
    values = pd.read_csv(GEN / fname).iloc[:, 0].dropna().astype(str).str.strip()
    out = set()
    for v in values:
        for part in re.split(r"\s+x\s+", v):  # intergenic "GENE1 x GENE2" -> both genes
            if part:
                out.add(harmonise(part))
    return out


def to_coding(x):
    if x in official:
        return x
    return prev.get(x)  # renamed protein-coding symbol, else None


G, C, S = {}, {}, {}
for d, (g, c) in FILES.items():
    G[d], C[d] = load(g), load(c)
    S[d] = G[d] | C[d]
cod = {d: {to_coding(x) for x in S[d] if to_coding(x)} for d in S}

# ---- counts per disease ----
rows = [dict(Disease=d, GWAS_Catalog=len(G[d]), ClinVar=len(C[d]) if FILES[d][1] else "NA",
             Both=len(G[d] & C[d]), Total_entities=len(S[d]), Protein_coding=len(cod[d])) for d in S]
counts = pd.DataFrame(rows)
counts.to_csv(OUT / "Table_genomics_counts.csv", index=False)
print(counts.to_string())


# ---- shared genes ----
def shared(sets, tag):
    allg = sorted(set().union(*sets.values()))
    M = pd.DataFrame({d: [g in sets[d] for g in allg] for d in sets}, index=allg)
    M["n"] = M.sum(axis=1)
    M["diseases"] = M.drop(columns="n").apply(lambda r: "+".join(r.index[r]), axis=1)
    sh = M[M.n >= 2].sort_values(["n", "diseases"], ascending=[False, True])
    sh.to_csv(OUT / f"shared_{tag}.csv")
    print(tag, "union", len(allg), M.n.value_counts().sort_index().to_dict())
    return M, sh


M_all, sh_all = shared(S, "all_entities")
M_cod, sh_cod = shared(cod, "coding_genes")

# ---- Venn diagrams ----
for sets, tag, title in [(S, "all_entities", "All gene entities"), (cod, "coding_genes", "Protein-coding genes")]:
    fig, ax = plt.subplots(figsize=(9, 8))
    venn({LAB[d]: s for d, s in sets.items()}, ax=ax, fontsize=10, legend_loc="upper left")
    ax.set_title(f"Overlap of genomic-layer {title.lower()} (GWAS Catalog + ClinVar)")
    fig.savefig(OUT / f"Venn_genomics_{tag}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"Venn_genomics_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)

# ---- permutation test (protein-coding genes; universe = HGNC protein-coding) ----
rng = np.random.default_rng(1)
universe = np.array(sorted(official))
sizes = [len(v) for v in cod.values()]
obs = {k: int((M_cod.n >= k).sum()) for k in (2, 3, 4, 5)}
N = 10000
sim = {k: [] for k in obs}
for _ in range(N):
    c = np.zeros(len(universe), int)
    for n in sizes:
        c[rng.choice(len(universe), n, replace=False)] += 1
    for k in obs:
        sim[k].append((c >= k).sum())
perm = {k: (float(np.mean(sim[k])), int(np.sum(np.array(sim[k]) >= obs[k]) + 1) / (N + 1)) for k in obs}
print("coding sizes", dict(zip(cod, sizes)), "observed", obs)
print("permutation (mean expected, P)", perm)

# ---- Excel workbook ----
with pd.ExcelWriter(OUT / "Genomics_layer_tables.xlsx") as w:
    counts.to_excel(w, sheet_name="Counts", index=False)
    sh_all.to_excel(w, sheet_name="Shared_all_entities")
    sh_cod.to_excel(w, sheet_name="Shared_coding")
    pd.DataFrame({d: pd.Series(sorted(S[d])) for d in S}).to_excel(w, sheet_name="Entity_lists_all", index=False)
    pd.DataFrame({d: pd.Series(sorted(cod[d])) for d in cod}).to_excel(w, sheet_name="Entity_lists_coding", index=False)

# ---- pairwise overlaps and genes shared by >=3 diseases ----
for a, b in itertools.combinations(cod, 2):
    print(a, b, "coding:", len(cod[a] & cod[b]), "all:", len(S[a] & S[b]))
print(sh_all[sh_all.n >= 3][["n", "diseases"]].to_string())
print(sh_cod[sh_cod.n >= 3][["n", "diseases"]].to_string())


def source(g, d):
    a, b = g in G[d], g in C[d]
    return "GWAS+ClinVar" if a and b else "GWAS" if a else "ClinVar" if b else ""


print("| Gene | n | " + " | ".join(FILES) + " |")
for g, r in sh_cod[sh_cod.n >= 3].iterrows():
    print("| " + g + " | " + str(r.n) + " | " + " | ".join(source(g, d) or "-" for d in FILES) + " |")
