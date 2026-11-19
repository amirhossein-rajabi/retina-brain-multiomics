#!/usr/bin/env bash
# Run the complete analysis from the repository root.
# Input files must be in place as described in data/README.md.
set -euo pipefail
cd "$(dirname "$0")"

Rscript scripts/00_install.R
Rscript scripts/01_metadata_explore.R GSE48350 GSE114517 GSE115828 GSE160306 GSE45570
Rscript scripts/02_DEG.R
Rscript scripts/03_shared_genes.R                 # FDR-based membership (sensitivity)
PCOL=pvalue Rscript scripts/03_shared_genes.R     # nominal P (primary analysis)
Rscript scripts/04_shared_neuro_ophth_2plus.R
python3 scripts/05_genomics_overlap.py
python3 scripts/06_genomics_and_supp.py
