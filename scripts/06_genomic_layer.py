"""Genomic layer: gene entities, shared genes, permutation test and summary tables.

Shared genes are present in >=2 diseases including at least one neurodegenerative
(AD, PD) and one ophthalmic (AMD, DR, glaucoma) disease. The script

* merges GWAS Catalog and ClinVar gene entities per disease and harmonises them
  to HGNC nomenclature (Table 1),
* identifies shared gene entities / protein-coding genes (Tables 2-3),
* runs a 10,000-iteration size-matched permutation test against the HGNC
  protein-coding universe,
* draws Venn diagrams of the genomic and transcriptomic overlaps, and
* writes an Excel workbook with the per-disease and shared gene lists.

Requires the outputs of scripts 02-04 in results/.

Usage (from the repository root):
    python scripts/06_genomic_layer.py
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
perm = {k: dict(observed=obs[k], expected_mean=round(float(np.mean(sim[k])), 2), permutation_max=int(np.max(sim[k])),
                P=(int(np.sum(np.array(sim[k]) >= obs[k])) + 1) / (N + 1)) for k in obs}


# ---- summary tables (Markdown) ----
def md_shared(M, gtitle):
    r = M[M.qual & (M.n >= 3)].sort_values(["n", "Diseases"], ascending=[False, True])
    return (f"| {gtitle} | Number of diseases | Diseases |\n|---|---|---|\n"
            + "\n".join(f"| {g} | {x.n} | {x.Diseases} |" for g, x in r.iterrows()))


t_all, t_cod = md_shared(MA, "Gene entity"), md_shared(MC, "Gene")

tx_shared_in = pd.read_excel(RES / "Shared_genes_neuro_ophthalmic_2plus.xlsx", sheet_name="Shared_genes")
t3 = tx_shared_in[tx_shared_in.n_diseases >= 3].copy()
t3["diseases"] = t3.diseases.str.replace("+", ", ").str.replace("GL", "Glaucoma")
t_tx = ("| Gene | Number of diseases | Diseases | Direction of change |\n|---|---|---|---|\n"
        + "\n".join(f"| {r.gene} | {r.n_diseases} | {r.diseases} | "
                    f"{'Same direction (up)' if r.direction == 'all up' else 'Discordant'} |" for r in t3.itertuples()))

count_rows = [(d, len(G[d]), len(C[d]) if FILES[d][1] else "NA", len(G[d] & C[d]) if FILES[d][1] else "NA",
               len(S[d]), len(cod[d])) for d in D]
t1 = ("| Disease | GWAS Catalog | ClinVar | Both sources | Total gene entities | Protein-coding genes |\n"
      "|---|---|---|---|---|---|\n"
      + "\n".join("| " + " | ".join(f"{x:,}" if isinstance(x, int) else str(x) for x in r) + " |" for r in count_rows))


# ---- Excel workbook ----
def entity_type(g):
    if g in official:
        return "Protein-coding gene"
    return "Intergenic region" if " x " in g else "Non-coding / other"


entities_df = pd.DataFrame(
    [(d, g, entity_type(g),
      "GWAS Catalog + ClinVar" if g in G[d] and g in C[d] else "GWAS Catalog" if g in G[d] else "ClinVar")
     for d in D for g in sorted(S[d])],
    columns=["Disease", "Gene entity", "Entity type", "Source"])


def shared_table(M):
    r = M[M.qual].sort_values(["n", "Diseases"], ascending=[False, True])
    r = r[["n", "Diseases"]].rename(columns={"n": "Number of diseases"})
    return r.reset_index(names="Gene entity")


shared_all_df, shared_coding_df = shared_table(MA), shared_table(MC)

tx_shared_df = tx_shared_in[["gene", "n_diseases", "diseases", "direction"]
         + [f"{d}_{m}" for d in TX_KEYS for m in ("log2FC", "pvalue")]].copy()
tx_shared_df["diseases"] = tx_shared_df.diseases.str.replace("+", ", ").str.replace("GL", "Glaucoma")
tx_shared_df.columns = [c.replace("GL_", "Glaucoma_").replace("_pvalue", " P").replace("_log2FC", " log2FC") for c in tx_shared_df.columns]
tx_shared_df = tx_shared_df.rename(columns={"gene": "Gene", "n_diseases": "Number of diseases", "diseases": "Diseases",
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
search_terms_df = pd.DataFrame([(d, *TERMS[d]) for d in D], columns=["Disease", "GWAS Catalog trait terms", "ClinVar condition terms"])

samples_df = pd.read_csv(RES / "TableS1_samples.csv").drop(columns=["included", "reason"])
samples_df["Disease"] = samples_df.Disease.replace({"Alzheimer's disease": "AD", "Parkinson's disease": "PD",
                                    "Diabetic retinopathy": "DR"})

sheets = {
    "samples": samples_df,                    # transcriptomic samples included
    "transcriptomic_shared": tx_shared_df,    # shared DEGs, direction and per-disease statistics
    "genomic_entities": entities_df,          # genomic entities per disease, entity type and source
    "genomic_shared_all": shared_all_df,      # shared genomic entities
    "genomic_shared_coding": shared_coding_df,  # shared protein-coding genes
    "search_terms": search_terms_df,          # GWAS Catalog / ClinVar search terms
}
with pd.ExcelWriter(OUT / "genomic_and_shared_gene_tables.xlsx") as w:
    for name, df in sheets.items():
        df.to_excel(w, sheet_name=name, index=False)


# ---- summary numbers ----
def n_exact(M, k):
    return int(((M.n == k) & M.qual).sum())


summary = dict(
    shared_all_entities=int(MA.qual.sum()), shared_coding=int(MC.qual.sum()),
    shared_by_n_all={k: n_exact(MA, k) for k in (2, 3, 4, 5)},
    shared_by_n_coding={k: n_exact(MC, k) for k in (2, 3, 4, 5)},
    permutation=perm, observed=obs, transcriptomic_DEGs_total_up_down=tx,
    union_all_entities=len(MA), union_coding=len(MC),
)
(OUT / "summary_numbers.json").write_text(json.dumps(summary, indent=2, default=str))
(OUT / "summary_tables.md").write_text("\n\n".join([t1, t_all, t_cod, t_tx]) + "\n")
print(json.dumps(summary, indent=2, default=str))
print(t1, t_all, t_cod, t_tx, sep="\n\n")
