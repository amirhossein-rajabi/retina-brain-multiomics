# Input data

## Genomic layer (included)

`data/genomics/` contains the gene lists used in the study, one file per disease and
source, with gene symbols in the first column. They were exported from the
[NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) (mapped genes of associations
with *P* < 5 × 10⁻⁸) and [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
(pathogenic / likely pathogenic, review status ≥ 2 stars), accessed on
12 November 2026. Both databases are updated continuously, so these files are kept
as a snapshot.

| Disease | GWAS Catalog | ClinVar |
|---|---|---|
| AD | `ADgenesGWAS.csv` | `ADgenesClinvar.csv` |
| PD | `PDgenesGWAS.csv` | `PDGenesClinvar.csv` |
| AMD | `AMDgenesGWAS.csv` | `ARMDgenesClinvar.csv` |
| DR | `DRgenesGWAS.csv` | no records met the criteria |
| Glaucoma | `GlaucomagenesGWAS.csv` | `GlaucomaGenesClinVar.csv` |

Intergenic entries are written as `GENE1 x GENE2`.

Search terms (*Homo sapiens*):

| Disease | GWAS Catalog traits | ClinVar conditions |
|---|---|---|
| AD | Alzheimer's disease; late-onset; early-onset | Alzheimer disease; familial Alzheimer disease |
| PD | Parkinson's disease; early-onset | Parkinson disease; early-onset Parkinson disease |
| AMD | age-related macular degeneration; neovascular (wet) AMD; geographic atrophy; advanced AMD | age-related macular degeneration |
| DR | diabetic retinopathy; proliferative DR; non-proliferative DR; diabetic macular edema | diabetic retinopathy; proliferative DR; non-proliferative DR |
| Glaucoma | glaucoma; primary open-angle glaucoma; normal-tension glaucoma | glaucoma; primary open-angle glaucoma; normal-tension glaucoma |

## Reference (included)

`data/reference/hgnc_protein_coding.txt` — snapshot of the HGNC protein-coding gene
set (19,276 genes; columns `hgnc_id`, `symbol`, `prev_symbol`) from
[genenames.org](https://www.genenames.org/download/statistics-and-files/). It is used
to harmonise symbols and as the universe of the genomic permutation test.

## Transcriptomic layer (downloaded)

The GEO files are not tracked. `scripts/01_download_geo.R` downloads them into
`data/`:

| Disease | Accession | Files |
|---|---|---|
| AD | [GSE48350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350) | series matrix |
| PD | [GSE114517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114517) | series matrix, `GSE114517_RAW.tar` → `data/PD/` |
| AMD | [GSE115828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115828) | series matrix, `GSE115828_RSEM_gene_counts.tsv.gz` → `data/AMD/` |
| DR | [GSE160306](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306) | series matrix, gene counts → `data/DR/GSE160306_counts.txt.gz` |
| Glaucoma | [GSE45570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45570) | series matrix |
