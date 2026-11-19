# Transcriptomic shared genes: >=2 diseases with >=1 neurodegenerative (AD/PD) and >=1 ophthalmic
# (AMD/DR/Glaucoma) disease; nominal P < 0.05 and |log2FC| > 0.5. Adds direction of change.
# Requires the nominal-P membership table from: PCOL=pvalue Rscript scripts/03_shared_genes.R
# Usage: Rscript scripts/04_shared_neuro_ophth_2plus.R
suppressMessages(library(writexl))
t <- read.csv("results/all_DEG_membership_LFC05_nominalP.csv")
d <- c("AD","PD","AMD","DR","GL")
s <- t[t$neuro & t$ophthalmic & t$n_diseases >= 2, ]
fc <- as.matrix(s[paste0(d, "_log2FC")]); hit <- as.matrix(s[d])
s$direction <- sapply(seq_len(nrow(s)), function(i) { x <- fc[i, hit[i, ]]; if (all(x > 0)) "all up" else if (all(x < 0)) "all down" else "discordant" })
s$criterion <- ifelse(s$n_diseases >= 3, ">=3 diseases", "2 diseases")
s <- s[order(-s$n_diseases, s$diseases, s$gene), ]
out <- s[c("gene","n_diseases","diseases","direction","criterion", as.vector(rbind(paste0(d,"_log2FC"), paste0(d,"_padj"))))]
pv <- function(dd) { r <- readRDS("results/DEG/all_DEG.rds")[[dd]]; r$pvalue[match(out$gene, r$gene)] }
for (dd in d) out[[paste0(dd, "_pvalue")]] <- signif(pv(dd), 3)
write_xlsx(list(Shared_genes = out,
  Notes = data.frame(Note = c("Significance: nominal P < 0.05 and |log2FC| > 0.5 in each dataset",
    "Inclusion: significant in >=2 diseases, at least one neurological (AD, PD) and one ophthalmic (AMD, DR, Glaucoma)",
    "AD=GSE48350 hippocampus; PD=GSE114517 substantia nigra; AMD=GSE115828 MGS4 vs MGS1; DR=GSE160306 macula NPDR+PDR vs control; GL=GSE45570 optic nerve head POAG vs control",
    "direction: concordant up/down across the diseases where the gene is significant; NA = gene not tested in that dataset"))),
  "results/Shared_genes_neuro_ophthalmic_2plus.xlsx")
cat("Total:", nrow(out), "| >=3 diseases:", sum(out$n_diseases >= 3), "| 2 diseases:", sum(out$n_diseases == 2), "\n")
print(table(out$diseases[out$n_diseases == 2], out$direction[out$n_diseases == 2]))
for (p in unique(out$diseases)) cat(p, ":", paste(out$gene[out$diseases == p], collapse = ", "), "\n")
