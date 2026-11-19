# Download the GEO data used in the transcriptomic layer.
# - series matrices and parsed sample metadata (data/<GSE>_meta.rds), with a
#   summary of the sample characteristics printed for each series
# - supplementary count files of the RNA-seq datasets (PD, AMD, DR)
# Usage: Rscript scripts/01_download_geo.R            # all five datasets
#        Rscript scripts/01_download_geo.R GSE48350   # selected series only
suppressMessages(library(GEOquery))
options(timeout = 3600)
dir.create("data", showWarnings = FALSE)

# Turn every "key: value" entry of the characteristics columns into its own column
parse_ch <- function(p) {
  cc <- grep("^characteristics", colnames(p), value = TRUE)
  out <- data.frame(gsm = rownames(p), title = p$title, source = p$source_name_ch1, platform = p$platform_id)
  for (c in cc) for (i in seq_len(nrow(p))) {
    v <- p[i, c]; if (!nzchar(v) || !grepl(":", v)) next
    k <- trimws(tolower(sub(":.*", "", v))); out[i, k] <- trimws(sub("^[^:]*:", "", v))
  }
  out
}

ids <- commandArgs(TRUE)
if (!length(ids)) ids <- c("GSE48350", "GSE114517", "GSE115828", "GSE160306", "GSE45570")
for (g in ids) {
  gs <- NULL
  for (attempt in 1:4) {
    gs <- tryCatch(suppressMessages(getGEO(g, GSEMatrix = TRUE, destdir = "data", getGPL = FALSE)),
                   error = function(e) NULL)
    if (!is.null(gs)) break
    Sys.sleep(5)
  }
  if (is.null(gs)) { cat(g, "FAILED\n"); next }
  m <- Reduce(function(a, b) merge(a, b, all = TRUE), lapply(gs, function(e) parse_ch(pData(e))))
  saveRDS(m, sprintf("data/%s_meta.rds", g))
  cat("\n=====", g, "| platforms", paste(sapply(gs, annotation), collapse = ","),
      "| samples", nrow(m), "| features", paste(sapply(gs, nrow), collapse = ","), "\n")
  for (k in setdiff(colnames(m), c("gsm", "title"))) {
    u <- table(m[[k]], useNA = "ifany")
    if (length(u) > 12) cat(" ", k, ": ", length(u), " unique, e.g. ", paste(head(names(u), 5), collapse = " | "), "\n", sep = "")
    else cat(" ", k, ": ", paste(names(u), u, sep = "=", collapse = " | "), "\n", sep = "")
  }
  cat("  example titles:", paste(head(m$title, 3), collapse = " | "), "\n")
}

# ---- supplementary count files (RNA-seq datasets) ----
geo_suppl <- function(gse, file) {
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/%s", substr(gse, 1, nchar(gse) - 3), gse, file)
}
fetch <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(dest)) download.file(url, dest, mode = "wb")
  cat("  ", dest, "\n")
}
cat("\n===== supplementary files\n")
if ("GSE114517" %in% ids) {                       # PD: per-sample count tables
  tar <- "data/PD/GSE114517_RAW.tar"
  fetch(geo_suppl("GSE114517", "GSE114517_RAW.tar"), tar)
  untar(tar, exdir = "data/PD")
}
if ("GSE115828" %in% ids)                         # AMD: RSEM gene counts
  fetch(geo_suppl("GSE115828", "GSE115828_RSEM_gene_counts.tsv.gz"), "data/AMD/GSE115828_RSEM_gene_counts.tsv.gz")
if ("GSE160306" %in% ids)                         # DR: gene counts
  fetch(geo_suppl("GSE160306", "GSE160306_human_retina_DR_totalRNA_counts.txt.gz"), "data/DR/GSE160306_counts.txt.gz")
