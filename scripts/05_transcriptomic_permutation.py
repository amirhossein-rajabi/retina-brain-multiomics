"""Permutation test for the transcriptomic shared genes.

Shared genes are DEGs (nominal P < 0.05, |log2FC| > 0.5) present in >=2 diseases,
including at least one neurodegenerative (AD, PD) and one ophthalmic (AMD, DR,
glaucoma) disease. In each of 10,000 iterations, a random gene set of the observed
DEG size is drawn without replacement from the genes tested in each dataset, and
the number of genes meeting the shared-gene definition in >=2 and >=3 diseases is
recorded. Empirical P = (r + 1) / (n + 1).

Requires results/DEG/<disease>_all_genes.csv from 02_DEG.R.

Usage (from the repository root):
    python scripts/05_transcriptomic_permutation.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEG_DIR = ROOT / "results" / "DEG"
OUT = ROOT / "results" / "transcriptomic_permutation.json"

DATASETS = {"AD": "AD", "PD": "PD", "AMD": "AMD", "DR": "DR", "Glaucoma": "Glaucoma"}
NEURO, OPHTH = ["AD", "PD"], ["AMD", "DR", "Glaucoma"]
P_CUT, LFC_CUT = 0.05, 0.5
N_PERM, SEED = 10000, 1

tested, degs = {}, {}
for disease, fname in DATASETS.items():
    res = pd.read_csv(DEG_DIR / f"{fname}_all_genes.csv").dropna(subset=["gene"])
    tested[disease] = sorted(set(res.gene))
    degs[disease] = set(res.gene[(res.pvalue < P_CUT) & (res.log2FC.abs() > LFC_CUT)])

keys = list(DATASETS)
genes = sorted(set().union(*tested.values()))
index = {g: i for i, g in enumerate(genes)}
neuro_cols = [keys.index(d) for d in NEURO]
ophth_cols = [keys.index(d) for d in OPHTH]


def shared_counts(B):
    """Number of genes meeting the shared-gene definition in >=2 and >=3 diseases."""
    n = B.sum(1)
    q = B[:, neuro_cols].any(1) & B[:, ophth_cols].any(1)
    return {k: int(((n >= k) & q).sum()) for k in (2, 3)}


B_obs = np.zeros((len(genes), len(keys)), bool)
for j, d in enumerate(keys):
    B_obs[[index[g] for g in degs[d]], j] = True
observed = shared_counts(B_obs)

rng = np.random.default_rng(SEED)
tested_idx = {d: np.array([index[g] for g in tested[d]]) for d in keys}
sim = {k: np.empty(N_PERM, int) for k in observed}
for i in range(N_PERM):
    B = np.zeros_like(B_obs)
    for j, d in enumerate(keys):
        B[rng.choice(tested_idx[d], len(degs[d]), replace=False), j] = True
    for k, v in shared_counts(B).items():
        sim[k][i] = v

summary = {
    "genes_tested": {d: len(tested[d]) for d in keys},
    "DEGs": {d: len(degs[d]) for d in keys},
    "iterations": N_PERM,
    "seed": SEED,
}
for k in observed:
    r = int((sim[k] >= observed[k]).sum())
    summary[f"shared_in_at_least_{k}"] = {
        "observed": observed[k],
        "expected_mean": round(float(sim[k].mean()), 2),
        "permutation_max": int(sim[k].max()),
        "P": (r + 1) / (N_PERM + 1),
    }

OUT.write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
