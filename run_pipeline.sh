#!/usr/bin/env bash
# Run the complete analysis from the repository root.
set -euo pipefail
cd "$(dirname "$0")"

Rscript scripts/00_install.R
Rscript scripts/01_download_geo.R
Rscript scripts/02_DEG.R
Rscript scripts/03_shared_genes.R                 # BH-FDR membership (sensitivity analysis)
PCOL=pvalue Rscript scripts/03_shared_genes.R     # nominal P membership (primary analysis)
Rscript scripts/04_shared_neuro_ophth_2plus.R
python3 scripts/05_transcriptomic_permutation.py
python3 scripts/06_genomic_layer.py
