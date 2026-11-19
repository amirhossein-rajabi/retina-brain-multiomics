# Sample selection, dataset characteristics (Table 4) and differential expression
# analysis for the five GEO datasets (limma for microarrays, DESeq2 for RNA-seq).
# Requires data/ prepared as described in data/README.md and 01_metadata_explore.R.
# Usage (from the repository root): Rscript scripts/02_DEG.R
suppressMessages({library(GEOquery); library(limma); library(DESeq2); library(AnnotationDbi)
  library(org.Hs.eg.db); library(hgu133plus2.db); library(huex10sttranscriptcluster.db); library(writexl)})
dir.create("results/DEG", recursive = TRUE, showWarnings = FALSE)
LFC <- as.numeric(Sys.getenv("LFC", "0.5")); FDR <- 0.05
meta_all <- list(); table1 <- list()

ens2sym <- function(ens) { ens <- sub("\\..*", "", ens)
  unname(mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")) }

# ---- generic analysis helpers ----
run_array <- function(expr, sym, grp, name) {
  if (max(expr, na.rm = TRUE) > 100) expr <- log2(pmax(expr, 1))
  expr <- normalizeBetweenArrays(expr, method = "quantile")
  ok <- !is.na(sym) & sym != ""; expr <- expr[ok, ]; sym <- sym[ok]
  keep <- tapply(seq_along(sym), sym, function(i) i[which.max(rowMeans(expr[i, , drop = FALSE]))])
  expr <- expr[unlist(keep), ]; rownames(expr) <- names(keep)
  d <- model.matrix(~ factor(grp, levels = c("Control", "Disease")))
  tt <- topTable(eBayes(lmFit(expr, d)), coef = 2, number = Inf)
  data.frame(gene = rownames(tt), log2FC = tt$logFC, pvalue = tt$P.Value, padj = tt$adj.P.Val)
}
run_seq <- function(counts, grp, name) {
  counts <- round(as.matrix(counts)); storage.mode(counts) <- "integer"
  keep <- rowSums(counts >= 10) >= min(table(grp)); counts <- counts[keep, ]
  dds <- DESeqDataSetFromMatrix(counts, data.frame(grp = factor(grp, levels = c("Control", "Disease"))), ~ grp)
  dds <- DESeq(dds, quiet = TRUE); r <- as.data.frame(results(dds, name = "grp_Disease_vs_Control"))
  r$gene <- ens2sym(rownames(r)); r <- r[!is.na(r$gene) & !is.na(r$padj), ]
  r <- r[order(r$pvalue), ]; r <- r[!duplicated(r$gene), ]
  data.frame(gene = r$gene, log2FC = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
}
summ_age <- function(x) { x <- as.numeric(x); x <- x[!is.na(x)]
  if (!length(x)) return("NR"); sprintf("%.1f ± %.1f (%g–%g)", mean(x), sd(x), min(x), max(x)) }
summ_sex <- function(s) { s <- toupper(substr(trimws(s), 1, 1)); if (all(is.na(s))) return("NR")
  sprintf("%d/%d", sum(s == "M", na.rm = TRUE), sum(s == "F", na.rm = TRUE)) }
add_table1 <- function(disease, gse, platform, tissue, comp, m, excl) {
  d <- m[m$group == "Disease", ]; c <- m[m$group == "Control", ]
  pa <- tryCatch(signif(wilcox.test(as.numeric(d$age), as.numeric(c$age), exact = FALSE)$p.value, 2), error = function(e) NA)
  ps <- tryCatch(signif(fisher.test(table(m$group, toupper(substr(m$sex, 1, 1))))$p.value, 2), error = function(e) NA)
  table1[[disease]] <<- data.frame(Disease = disease, GEO = gse, Platform = platform, Tissue = tissue,
    Comparison = comp, n_Disease = nrow(d), n_Control = nrow(c), Excluded = excl,
    Age_Disease = summ_age(d$age), Age_Control = summ_age(c$age), p_age = pa,
    Sex_MF_Disease = summ_sex(d$sex), Sex_MF_Control = summ_sex(c$sex), p_sex = ps, check.names = FALSE)
  meta_all[[disease]] <<- cbind(Disease = disease, m)
}
save_deg <- function(res, name) {
  write.csv(res, sprintf("results/DEG/%s_all_genes.csv", name), row.names = FALSE)
  sig <- res[res$padj < FDR & abs(res$log2FC) > LFC, ]
  cat(sprintf("%-10s genes tested=%5d  DEG=%4d (up %d / down %d)\n", name, nrow(res), nrow(sig),
              sum(sig$log2FC > 0), sum(sig$log2FC < 0)))
  res
}
deg <- list()

# ---- 1. Alzheimer: GSE48350, hippocampus, AD vs controls aged >=60 ----
e <- getGEO(filename = "data/GSE48350_series_matrix.txt.gz", getGPL = FALSE)
m <- readRDS("data/GSE48350_meta.rds"); m <- m[m[["brain region"]] == "hippocampus", ]
m$age <- as.numeric(m[["age (yrs)"]]); m$sex <- m$gender
m$group <- ifelse(!is.na(m[["braak stage"]]), "Disease", "Control")
n_young <- sum(m$group == "Control" & m$age < 60)
m$included <- !(m$group == "Control" & m$age < 60); m$reason <- ifelse(m$included, "", "Control <60 y")
ms <- m[m$included, ]
add_table1("Alzheimer's disease", "GSE48350", "Affymetrix HG-U133 Plus 2.0 (GPL570)", "Hippocampus",
           "AD vs non-demented controls", ms, sprintf("%d controls aged <60 y; other brain regions", n_young))
x <- exprs(e)[, ms$gsm]
deg$AD <- save_deg(run_array(x, mapIds(hgu133plus2.db, rownames(x), "SYMBOL", "PROBEID"), ms$group), "AD")

# ---- 2. Parkinson: GSE114517, substantia nigra, PD vs control ----
m <- readRDS("data/GSE114517_meta.rds"); m <- m[m$tissue == "Substantia nigra", ]
m$group <- ifelse(grepl("PD", m$title), "Disease", "Control")
m$sex <- ifelse(grepl("female", m$title, ignore.case = TRUE), "F", "M")
m$included <- TRUE; m$reason <- ""
add_table1("Parkinson's disease", "GSE114517", "RNA-seq, Illumina NextSeq 500 (GPL18573)", "Substantia nigra",
           "PD vs non-demented controls", m, "Amygdala and medial temporal gyrus samples")
cnt <- sapply(m$gsm, function(g) { f <- list.files("data/PD", paste0("^", g, "_.*tabular"), full.names = TRUE)
  t <- read.delim(f, header = FALSE); setNames(t[[2]], t[[1]]) })
cnt <- cnt[grepl("^ENSG", rownames(cnt)), ]
deg$PD <- save_deg(run_seq(cnt, m$group), "PD")

# ---- 3. AMD: GSE115828, MGS4 vs MGS1 ----
m <- readRDS("data/GSE115828_meta.rds"); m <- m[m$mgs_level %in% c("1", "4"), ]
m$group <- ifelse(m$mgs_level == "4", "Disease", "Control"); m$included <- TRUE; m$reason <- ""
add_table1("AMD", "GSE115828", "RNA-seq, Illumina HiSeq 2500 (GPL16791)", "Retina (macular punch)",
           "MGS4 vs MGS1", m, "MGS2, MGS3 and ungraded samples")
cnt <- read.delim("data/AMD/GSE115828_RSEM_gene_counts.tsv.gz", row.names = 1, check.names = FALSE)
stopifnot(all(m$title %in% colnames(cnt)))
deg$AMD <- save_deg(run_seq(cnt[, m$title], m$group), "AMD")

# ---- 4. Diabetic retinopathy: GSE160306, macula, NPDR+PDR vs healthy controls ----
m <- readRDS("data/GSE160306_meta.rds"); m <- m[m$sample_site == "Macula" & m$disease_group != "Diabetic", ]
m$group <- ifelse(m$disease_group == "Control", "Control", "Disease"); m$sex <- m$gender
flag <- grepl("glaucoma|AMD|ARMD|drusen", paste(m$pathology_note, m$note_other), ignore.case = TRUE)
m$included <- !flag; m$reason <- ifelse(flag, paste("AMD/glaucoma suspect:", m$pathology_note, m$note_other), "")
ms <- m[m$included, ]
cnt <- read.csv("data/DR/GSE160306_counts.txt.gz", row.names = 1, check.names = FALSE)[, ms$sampleid]
cnt <- t(rowsum(t(cnt), ms$donor))                       # merge OD/OS eyes of the same donor
md <- ms[!duplicated(ms$donor), ]; cnt <- cnt[, md$donor]
add_table1("Diabetic retinopathy", "GSE160306", "RNA-seq, Illumina HiSeq 4000 (GPL20301)", "Retina (macula)",
           "NPDR + PDR vs healthy controls", md,
           sprintf("Diabetics without DR; peripheral retina; %d macular samples with AMD/glaucoma suspect", sum(flag)))
deg$DR <- save_deg(run_seq(cnt, md$group), "DR")

# ---- 5. Glaucoma: GSE45570, POAG vs unaffected control ----
e <- getGEO(filename = "data/GSE45570_series_matrix.txt.gz", getGPL = FALSE)
m <- readRDS("data/GSE45570_meta.rds"); m <- m[m$group != "Ocular Hypertension", ]
m$group <- ifelse(m$group == "Glaucoma Affected", "Disease", "Control")
m$age <- as.numeric(sub(" yrs old", "", m$age)); m$sex <- NA; m$included <- TRUE; m$reason <- ""
add_table1("Glaucoma", "GSE45570", "Affymetrix Human Exon 1.0 ST (GPL5175)", "Optic nerve head",
           "POAG vs unaffected controls", m, "Ocular hypertension samples")
x <- exprs(e)[, m$gsm]
deg$GL <- save_deg(run_array(x, mapIds(huex10sttranscriptcluster.db, rownames(x), "SYMBOL", "PROBEID"), m$group), "Glaucoma")

# ---- save ----
t1 <- do.call(rbind, table1); rownames(t1) <- NULL
write.csv(t1, "results/Table1_metadata.csv", row.names = FALSE)
keepc <- c("Disease", "gsm", "title", "group", "age", "sex", "included", "reason")
s1 <- do.call(rbind, lapply(meta_all, function(d) { for (k in keepc) if (!k %in% names(d)) d[[k]] <- NA; d[, keepc] }))
write.csv(s1, "results/TableS1_samples.csv", row.names = FALSE)
saveRDS(deg, "results/DEG/all_DEG.rds")
print(t1[, c(1, 6, 7, 9, 10, 11, 12, 13, 14)], row.names = FALSE)
