args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE))))
root <- dirname(here)
take <- list(
  "scripts/12_run_arms_k48all.R" = c("read_bin", "add_metrics", "drop_truncated", "run_mad",
                                     "run_spotsweeper", "check_source", "read_full_array",
                                     "aggregate_full_array", "run_spatialartifacts_v", "project"),
  "R/ad.R" = c("barcode_based_deconvolve"))
out <- character(0)
for (src in names(take)) {
  ex <- parse(file.path(root, src), keep.source = TRUE)
  sr <- attr(ex, "srcref")
  nm <- vapply(ex, function(e) if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
                                     is.name(e[[2]])) as.character(e[[2]]) else "", "")
  for (f in take[[src]]) {
    i <- which(nm == f)
    if (length(i) != 1) stop(sprintf("%s: found %d definitions of %s", src, length(i), f))
    out <- c(out, as.character(sr[[i]]), "")
  }
}
writeLines(out, file.path(here, "qc_core.R"))
cat("wrote", file.path(here, "qc_core.R"), "\n")
