"""Genomic-layer shared genes, permutation test, figures and supplementary tables.

Shared genes are present in >=2 diseases including at least one neurodegenerative
(AD, PD) and one ophthalmic (AMD, DR, glaucoma) disease. The script

* merges GWAS Catalog and ClinVar gene entities per disease and harmonises them
  to HGNC nomenclature (Table 1),
* identifies shared gene entities / protein-coding genes (Tables 2-3),
* runs a 10,000-iteration size-matched permutation test against the HGNC
  protein-coding universe,
* draws the genomic and transcriptomic Venn diagrams, and
* writes the supplementary workbook (Tables S1-S5).

Requires the outputs of scripts 02-04 in results/.

Usage (from the repository root):
    python scripts/06_genomics_and_supp.py
"""
import json
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
RES = ROOT / "results"                       # transcriptomic results (scripts 02-04)
OUT = ROOT / "results" / "genomics"
(OUT / "figs").mkdir(parents=True, exist_ok=True)

D = ["AD", "PD", "AMD", "DR", "Glaucoma"]
NEU, OPH = {"AD", "PD"}, {"AMD", "DR", "Glaucoma"}
TX_KEYS = ["AD", "PD", "AMD", "DR", "GL"]   # column keys used by the R scripts
FILES = {
    "AD": ("ADgenesGWAS.csv", "ADgenesClinvar.csv"),
    "PD": ("PDgenesGWAS.csv", "PDGenesClinvar.csv"),
    "AMD": ("AMDgenesGWAS.csv", "ARMDgenesClinvar.csv"),
    "DR": ("DRgenesGWAS.csv", None),
    "Glaucoma": ("GlaucomagenesGWAS.csv", "GlaucomaGenesClinVar.csv"),
}
LAB = {"AD": "Alzheimer's disease", "PD": "Parkinson's disease", "AMD": "AMD",
       "DR": "Diabetic retinopathy", "Glaucoma": "Glaucoma"}

# ---- HGNC harmonisation ----
hg = pd.read_csv(HGNC, sep="\t", dtype=str)
official = set(hg.symbol)
prev = {}
for sym, p in zip(hg.symbol, hg.prev_symbol.fillna("")):
    for q in p.split("|"):
        if q:
            prev[q] = sym


def harmonise(x):
    return prev[x] if (x not in official and x in prev) else x


def load(fname):
    """Read a gene list; intergenic entries ("GENE1 x GENE2") are kept as single entities."""
    if fname is None:
        return set()
    values = pd.read_csv(GEN / fname).iloc[:, 0].dropna().astype(str).str.strip()
    return {v if " x " in v else harmonise(v) for v in values if v}


G = {d: load(FILES[d][0]) for d in D}
C = {d: load(FILES[d][1]) for d in D}
S = {d: G[d] | C[d] for d in D}
cod = {d: {x for x in S[d] if x in official} for d in D}


def membership(sets):
    allg = sorted(set().union(*sets.values()))
    M = pd.DataFrame({d: [g in sets[d] for g in allg] for d in D}, index=allg)
    M["n"] = M[D].sum(axis=1)
    M["Diseases"] = M[D].apply(lambda r: ", ".join(r.index[r.values]), axis=1)
    M["qual"] = (M.n >= 2) & M[list(NEU)].any(axis=1) & M[list(OPH)].any(axis=1)
    return M


MA, MC = membership(S), membership(cod)


def venn_plot(sets, title, name):
    fig, ax = plt.subplots(figsize=(9, 8))
    venn({LAB[d]: sets[d] for d in D}, ax=ax, fontsize=10, legend_loc="upper left")
    ax.set_title(title)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / "figs" / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


venn_plot(S, "Overlap of genomic-layer gene entities (GWAS Catalog + ClinVar)", "Fig_Venn_genomic_all_entities")
venn_plot(cod, "Overlap of genomic-layer protein-coding genes (GWAS Catalog + ClinVar)", "Fig_Venn_genomic_coding")

# ---- transcriptomic layer (from 03_shared_genes.R, nominal P) ----
tm = pd.read_csv(RES / "all_DEG_membership_LFC05_nominalP.csv")
tsets = {d: set(tm.gene[tm[k]]) for d, k in zip(D, TX_KEYS)}
venn_plot(tsets, "Overlap of differentially expressed genes (nominal P < 0.05, |log2FC| > 0.5)",
          "Fig_Venn_transcriptomic_DEGs")
tx = {}
for d, k in zip(D, TX_KEYS):
    fc = tm[tm[k]][k + "_log2FC"]
    tx[d] = (len(fc), int((fc > 0).sum()), int((fc < 0).sum()))

# ---- permutation test (criterion-based, HGNC protein-coding universe) ----
rng = np.random.default_rng(1)
universe = sorted(official)
sizes = [len(cod[d]) for d in D]
N = 10000
ni = [D.index(d) for d in D if d in NEU]
oi = [D.index(d) for d in D if d in OPH]
obs = {k: int(((MC.n >= k) & MC.qual).sum()) for k in (2, 3)}
sim = {k: [] for k in obs}
for _ in range(N):
    B = np.zeros((len(universe), len(D)), bool)
    for j, n in enumerate(sizes):
        B[rng.choice(len(universe), n, replace=False), j] = True
    n_ = B.sum(1)
    q = B[:, ni].any(1) & B[:, oi].any(1)
    for k in obs:
        sim[k].append(int(((n_ >= k) & q).sum()))
perm = {k: (float(np.mean(sim[k])), (int(np.sum(np.array(sim[k]) >= obs[k])) + 1) / (N + 1)) for k in obs}


# ---- manuscript tables (Markdown) ----
def md_shared(M, gtitle):
    r = M[M.qual & (M.n >= 3)].sort_values(["n", "Diseases"], ascending=[False, True])
    return (f"| {gtitle} | Number of diseases | Diseases |\n|---|---|---|\n"
            + "\n".join(f"| {g} | {x.n} | {x.Diseases} |" for g, x in r.iterrows()))


t_all, t_cod = md_shared(MA, "Gene entity"), md_shared(MC, "Gene")

sup = pd.read_excel(RES / "Shared_genes_neuro_ophthalmic_2plus.xlsx", sheet_name="Shared_genes")
t3 = sup[sup.n_diseases >= 3].copy()
t3["diseases"] = t3.diseases.str.replace("+", ", ").str.replace("GL", "Glaucoma")
t_tx = ("| Gene | Number of diseases | Diseases | Direction of change |\n|---|---|---|---|\n"
        + "\n".join(f"| {r.gene} | {r.n_diseases} | {r.diseases} | "
                    f"{'Same direction (up)' if r.direction == 'all up' else 'Discordant'} |" for r in t3.itertuples()))

count_rows = [(d, len(G[d]), len(C[d]) if FILES[d][1] else "NA", len(G[d] & C[d]) if FILES[d][1] else "NA",
               len(S[d]), len(cod[d])) for d in D]
t1 = ("| Disease | GWAS Catalog | ClinVar | Both sources | Total gene entities | Protein-coding genes |\n"
      "|---|---|---|---|---|---|\n"
      + "\n".join("| " + " | ".join(f"{x:,}" if isinstance(x, int) else str(x) for x in r) + " |" for r in count_rows))


# ---- supplementary workbook ----
def entity_type(g):
    if g in official:
        return "Protein-coding gene"
    return "Intergenic region" if " x " in g else "Non-coding / other"


S3 = pd.DataFrame(
    [(d, g, entity_type(g),
      "GWAS Catalog + ClinVar" if g in G[d] and g in C[d] else "GWAS Catalog" if g in G[d] else "ClinVar")
     for d in D for g in sorted(S[d])],
    columns=["Disease", "Gene entity", "Entity type", "Source"])


def shared_table(M):
    r = M[M.qual].sort_values(["n", "Diseases"], ascending=[False, True])
    r = r[["n", "Diseases"]].rename(columns={"n": "Number of diseases"})
    return r.reset_index(names="Gene entity")


S4a, S4b = shared_table(MA), shared_table(MC)

S2 = sup[["gene", "n_diseases", "diseases", "direction"]
         + [f"{d}_{m}" for d in TX_KEYS for m in ("log2FC", "pvalue")]].copy()
S2["diseases"] = S2.diseases.str.replace("+", ", ").str.replace("GL", "Glaucoma")
S2.columns = [c.replace("GL_", "Glaucoma_").replace("_pvalue", " P").replace("_log2FC", " log2FC") for c in S2.columns]
S2 = S2.rename(columns={"gene": "Gene", "n_diseases": "Number of diseases", "diseases": "Diseases",
                        "direction": "Direction"})

TERMS = {
    "AD": ("Alzheimer's disease; late-onset Alzheimer's disease; early-onset Alzheimer's disease",
           "Alzheimer disease; familial Alzheimer disease"),
    "PD": ("Parkinson's disease; early-onset Parkinson's disease",
           "Parkinson disease; early-onset Parkinson disease"),
    "AMD": ("age-related macular degeneration; neovascular (wet) AMD; geographic atrophy; advanced AMD",
            "age-related macular degeneration; macular degeneration, age-related"),
    "DR": ("diabetic retinopathy; proliferative diabetic retinopathy; non-proliferative diabetic retinopathy; "
           "diabetic macular edema",
           "diabetic retinopathy; proliferative diabetic retinopathy; non-proliferative diabetic retinopathy "
           "(no records met the inclusion criteria)"),
    "Glaucoma": ("glaucoma; primary open-angle glaucoma; normal-tension glaucoma",
                 "glaucoma; primary open-angle glaucoma; normal-tension glaucoma"),
}
S5 = pd.DataFrame([(d, *TERMS[d]) for d in D], columns=["Disease", "GWAS Catalog trait terms", "ClinVar condition terms"])

S1 = pd.read_csv(RES / "TableS1_samples.csv").drop(columns=["included", "reason"])
S1["Disease"] = S1.Disease.replace({"Alzheimer's disease": "AD", "Parkinson's disease": "PD",
                                    "Diabetic retinopathy": "DR"})

index = pd.DataFrame({
    "Table": ["S1", "S2", "S3", "S4", "S5"],
    "Content": [
        "Transcriptomic samples included in the analysis",
        "Genes shared by at least two diseases in the transcriptomic layer (one neurodegenerative, one ophthalmic)",
        "Genomic-layer gene entities per disease with entity type and source",
        "Genomic-layer gene entities shared by at least two diseases (one neurodegenerative, one ophthalmic); "
        "sheet 'S4a' all entities, 'S4b' protein-coding genes",
        "Search terms used for the GWAS Catalog and ClinVar",
    ]})
with pd.ExcelWriter(OUT / "Supplementary_Tables.xlsx") as w:
    index.to_excel(w, sheet_name="Index", index=False)
    for name, df in [("S1", S1), ("S2", S2), ("S3", S3), ("S4a", S4a), ("S4b", S4b), ("S5", S5)]:
        df.to_excel(w, sheet_name=name, index=False)


# ---- summary numbers ----
def n_exact(M, k):
    return int(((M.n == k) & M.qual).sum())


summary = dict(
    shared_all_entities=int(MA.qual.sum()), shared_coding=int(MC.qual.sum()),
    shared_by_n_all={k: n_exact(MA, k) for k in (2, 3, 4, 5)},
    shared_by_n_coding={k: n_exact(MC, k) for k in (2, 3, 4, 5)},
    permutation_expected_and_P=perm, observed=obs, transcriptomic_DEGs_total_up_down=tx,
    union_all_entities=len(MA), union_coding=len(MC),
)
(OUT / "summary_numbers.json").write_text(json.dumps(summary, indent=2, default=str))
(OUT / "manuscript_tables.md").write_text("\n\n".join([t1, t_all, t_cod, t_tx]) + "\n")
print(json.dumps(summary, indent=2, default=str))
print(t1, t_all, t_cod, t_tx, sep="\n\n")
