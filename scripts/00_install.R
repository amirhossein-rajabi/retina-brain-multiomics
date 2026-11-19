# Install the R/Bioconductor packages used by the R scripts.
# Usage: Rscript scripts/00_install.R
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("ggplot2", "ggVennDiagram", "writexl")
bioc <- c("GEOquery", "limma", "DESeq2", "AnnotationDbi", "org.Hs.eg.db",
          "hgu133plus2.db", "huex10sttranscriptcluster.db")
need_c <- cran[!sapply(cran, requireNamespace, quietly = TRUE)]
if (length(need_c)) install.packages(need_c)
need_b <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(need_b)) BiocManager::install(need_b, ask = FALSE, update = FALSE)
for (p in c(cran, bioc)) cat(sprintf("%-30s %s\n", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING"))
