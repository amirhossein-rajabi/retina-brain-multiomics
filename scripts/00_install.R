# Install the CRAN and Bioconductor packages required by the R scripts.
# Usage: Rscript scripts/00_install.R
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("data.table", "ggplot2", "VennDiagram", "ggVennDiagram", "UpSetR", "writexl", "R.utils")
bioc <- c("GEOquery", "limma", "edgeR", "DESeq2", "AnnotationDbi", "org.Hs.eg.db",
          "hgu133plus2.db", "huex10sttranscriptcluster.db")
need_c <- cran[!sapply(cran, requireNamespace, quietly = TRUE)]
if (length(need_c)) install.packages(need_c, type = "binary")
need_b <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(need_b)) BiocManager::install(need_b, ask = FALSE, update = FALSE, type = "binary")
for (p in c(cran, bioc)) cat(sprintf("%-30s %s\n", p, requireNamespace(p, quietly = TRUE)))
