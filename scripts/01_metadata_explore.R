# Download GEO series matrices, parse sample characteristics and summarise metadata.
# Saves data/<GSE>_meta.rds for each accession (used by 02_DEG.R).
# Usage: Rscript scripts/01_metadata_explore.R GSE48350 GSE114517 GSE115828 GSE160306 GSE45570
suppressMessages(library(GEOquery))
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
