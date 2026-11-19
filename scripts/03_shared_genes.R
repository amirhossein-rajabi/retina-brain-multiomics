# DEG membership across diseases and shared DEGs (>=3 diseases, with >=1 neurodegenerative
# (AD/PD) and >=1 ophthalmic (AMD/DR/Glaucoma) disease), plus a Venn diagram.
# Usage:
#   Rscript scripts/03_shared_genes.R                  # FDR < 0.05
#   PCOL=pvalue Rscript scripts/03_shared_genes.R      # nominal P < 0.05 (used in the manuscript)
#   LFC=0.5 sets the |log2FC| threshold (default 0.5)
suppressMessages({library(ggVennDiagram); library(ggplot2); library(writexl)})
LFC <- as.numeric(Sys.getenv("LFC", "0.5")); FDR <- 0.05; PCOL <- Sys.getenv("PCOL", "padj")  # "padj" (FDR) or "pvalue" (nominal)
deg <- readRDS("results/DEG/all_DEG.rds")
lab <- c(AD = "Alzheimer", PD = "Parkinson", AMD = "AMD", DR = "Diabetic retinopathy", GL = "Glaucoma")
sets <- lapply(deg, function(r) r$gene[r[[PCOL]] < FDR & abs(r$log2FC) > LFC])
cat("LFC >", LFC, "| DEGs:", paste(names(sets), lengths(sets), collapse = ", "), "\n")

genes <- unique(unlist(sets))
tab <- data.frame(gene = genes)
for (d in names(deg)) {
  r <- deg[[d]][match(genes, deg[[d]]$gene), ]
  tab[[paste0(d, "_log2FC")]] <- signif(r$log2FC, 3); tab[[paste0(d, "_padj")]] <- signif(r$padj, 3)
  tab[[d]] <- genes %in% sets[[d]]
}
tab$n_diseases <- rowSums(tab[names(deg)])
tab$diseases <- apply(tab[names(deg)], 1, function(x) paste(names(deg)[x], collapse = "+"))
tab$neuro <- tab$AD | tab$PD; tab$ophthalmic <- tab$AMD | tab$DR | tab$GL
dir_of <- function(i) { fc <- unlist(tab[i, paste0(names(deg)[unlist(tab[i, names(deg)])], "_log2FC")])
  if (all(fc > 0)) "all up" else if (all(fc < 0)) "all down" else "mixed" }
shared <- tab[tab$n_diseases >= 3 & tab$neuro & tab$ophthalmic, ]
if (nrow(shared)) { shared$direction <- sapply(as.integer(rownames(shared)), dir_of)
  shared <- shared[order(-shared$n_diseases, shared$gene), ] }
tag <- paste0(sub("\\.", "", sprintf("LFC%s", LFC)), ifelse(PCOL == "padj", "_FDR", "_nominalP"))
write.csv(tab, sprintf("results/all_DEG_membership_%s.csv", tag), row.names = FALSE)
write.csv(shared, sprintf("results/shared_genes_%s.csv", tag), row.names = FALSE)
cat("Shared genes (>=3 diseases, neuro+ophthalmic):", nrow(shared), "\n")
if (nrow(shared)) print(shared[, c("gene", "diseases", "direction")], row.names = FALSE)

names(sets) <- lab[names(sets)]
p <- ggVennDiagram(sets, label = "count", label_alpha = 0, set_size = 4) +
  scale_fill_gradient(low = "#F4F6FA", high = "#3A6EA5") + theme(legend.position = "none") +
  labs(title = sprintf("Overlap of DEGs (%s < 0.05, |log2FC| > %s)", ifelse(PCOL == "padj", "FDR", "P"), LFC)) +
  scale_x_continuous(expand = expansion(mult = 0.15))
ggsave(sprintf("results/Venn_%s.png", tag), p, width = 9, height = 8, dpi = 300)
ggsave(sprintf("results/Venn_%s.pdf", tag), p, width = 9, height = 8)
