# Input data

Raw data are not tracked in this repository. Place the files below in `data/`
(paths are relative to the repository root) before running the pipeline.

## Transcriptomic layer (GEO)

Series matrices and the parsed metadata (`<GSE>_meta.rds`) are created by
`scripts/01_metadata_explore.R`. Supplementary count files must be downloaded
from the GEO pages of each series.

| Disease | Accession | Files expected in `data/` | Source |
|---|---|---|---|
| AD | [GSE48350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350) | `GSE48350_series_matrix.txt.gz`, `GSE48350_meta.rds` | `01_metadata_explore.R` |
| PD | [GSE114517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114517) | `GSE114517_meta.rds`, `PD/GSM*_*.tabular.txt.gz` | `GSE114517_RAW.tar` (extract into `data/PD/`) |
| AMD | [GSE115828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115828) | `GSE115828_meta.rds`, `AMD/GSE115828_RSEM_gene_counts.tsv.gz` | GEO supplementary file |
| DR | [GSE160306](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306) | `GSE160306_meta.rds`, `DR/GSE160306_counts.txt.gz` | GEO supplementary file |
| Glaucoma | [GSE45570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45570) | `GSE45570_series_matrix.txt.gz`, `GSE45570_meta.rds` | `01_metadata_explore.R` |

## Genomic layer (GWAS Catalog and ClinVar)

One gene list per disease and source, with gene symbols in the first column.
Lists were exported from the [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/)
(mapped genes of associations with P < 5 × 10⁻⁸) and
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (pathogenic / likely pathogenic,
review status ≥ 2 stars). Search terms are listed in Supplementary Table S5.

```
data/genomics/
├── ADgenesGWAS.csv          ADgenesClinvar.csv
├── PDgenesGWAS.csv          PDGenesClinvar.csv
├── AMDgenesGWAS.csv         ARMDgenesClinvar.csv
├── DRgenesGWAS.csv          (no ClinVar records)
└── GlaucomagenesGWAS.csv    GlaucomaGenesClinVar.csv
```

Intergenic entries are written as `GENE1 x GENE2`.

## Reference

`data/reference/hgnc_protein_coding.txt` — HGNC protein-coding gene set
(tab-separated, with the columns `symbol` and `prev_symbol`), available from
[genenames.org](https://www.genenames.org/download/statistics-and-files/)
("Locus groups → protein-coding genes").

The manuscript used the HGNC release accessed on 14 November 2026 (19,276
protein-coding genes). HGNC is updated continuously, so a later release can shift
counts slightly (for example, a locus newly classified as protein-coding).
