suppressPackageStartupMessages({
  library(Matrix); library(arrow); library(jsonlite); library(SpatialExperiment)
  library(scuttle); library(SpatialArtifacts); library(visiumHDResDeconv)
})
args <- commandArgs(trailingOnly = TRUE)
cell <- args[1]; bin <- as.integer(args[2]); outdir <- args[3]
conds <- strsplit(args[4], ",")[[1]]
verify <- length(args) >= 5 && args[5] == "verify"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ex <- parse("/work/scripts/12_run_arms_k48all.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) eval(e, globalenv())
}
SOURCE_8UM <- Sys.getenv("SOURCE_8UM", "/src/square_008um")
CACHE <- Sys.getenv("SA_CACHE", "/cache")
READERS <- as.integer(Sys.getenv("SA_READERS", "2"))
say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }

spatial_steps <- function(spe, outlier, resolution) {
  bin_size_um <- if (resolution == "8um") 8 else 16
  buffer_width_bins <- round(80 / bin_size_um)
  min_cluster_size_bins <- ceiling(1280 / bin_size_um^2)
  x <- spe$array_col; y <- spe$array_row
  mnx <- min(x); mxx <- max(x); mny <- min(y); mxy <- max(y)
  in_buf <- function(x, y) x <= mnx + buffer_width_bins | x >= mxx - buffer_width_bins |
    y <= mny + buffer_width_bins | y >= mxy - buffer_width_bins
  spe$edge_artifact_edge <- in_buf(x, y) & outlier
  spe$edge_artifact_problem_id <- NA_character_
  spe$edge_artifact_problem_size <- 0
  xyz <- data.frame(x = x, y = y, outlier = outlier, stringsAsFactors = FALSE)
  rownames(xyz) <- colnames(spe)
  cl <- SpatialArtifacts:::problemAreas_WithMorphology_terra(.xyz = xyz, uniqueIdentifier = "s",
                                                              min_cluster_size = min_cluster_size_bins,
                                                              resolution = resolution)
  if (nrow(cl) > 0) {
    interior <- character(0)
    for (cid in unique(cl$clumpID)) {
      sp <- cl$spotcode[cl$clumpID == cid]
      j <- match(sp, colnames(spe))
      if (!in_buf(median(x[j]), median(y[j]))) interior <- c(interior, cid)
    }
    ia <- cl[cl$clumpID %in% interior, , drop = FALSE]
    if (nrow(ia) > 0) {
      j <- match(ia$spotcode, colnames(spe))
      spe$edge_artifact_problem_id[j] <- ia$clumpID
      spe$edge_artifact_problem_size[j] <- ia$clumpSize
    }
  }
  cls <- classifyEdgeArtifacts(spe, min_spots = 20, name = "edge_artifact")
  v <- colData(cls)[["edge_artifact_classification"]]
  cls$barcode[which(!is.na(v) & v != "not_artifact" & as.logical(cls$in_tissue))]
}

todo <- if (verify) conds else conds[!file.exists(file.path(outdir, paste0(conds, ".txt")))]
if (length(todo) == 0) quit(save = "no")
cd0 <- file.path(cell, todo[1])
check_source(cd0)

read_slot <- function(n = READERS) {
  repeat {
    for (k in seq_len(n)) {
      l <- filelock::lock(file.path(CACHE, sprintf("sa_read_%d.lock", k)), timeout = 0)
      if (!is.null(l)) return(l)
    }
    Sys.sleep(20)
  }
}
key <- file.path(CACHE, sprintf("%s_%dum", basename(cell), bin))
dir.create(CACHE, showWarnings = FALSE)
cache <- paste0(key, ".rds")
cell_lock <- filelock::lock(paste0(key, ".lock"))
if (file.exists(cache)) {
  z <- readRDS(cache); base <- z$base; skel <- z$skel; rm(z)
} else {
  slot <- read_slot()
  ft_d <- read.delim(gzfile(file.path(cd0, "square_008um", "filtered_feature_bc_matrix", "features.tsv.gz")), header = FALSE)
  mt <- grepl("^MT-", ft_d[[2]], ignore.case = TRUE)
  raw <- file.path(SOURCE_8UM, "raw_feature_bc_matrix")
  m_all <- Matrix::readMM(gzfile(file.path(raw, "matrix.mtx.gz")))
  bc_all <- readLines(gzfile(file.path(raw, "barcodes.tsv.gz")))
  ft <- read.delim(gzfile(file.path(raw, "features.tsv.gz")), header = FALSE)
  gi <- match(ft_d[[1]], ft[[1]])
  stopifnot(!anyNA(gi))
  m_all <- as(m_all[gi, , drop = FALSE], "CsparseMatrix")
  pos <- as.data.frame(arrow::read_parquet(file.path(SOURCE_8UM, "spatial", "tissue_positions.parquet")))
  pos <- pos[match(bc_all, pos$barcode), ]
  spe8 <- SpatialExperiment(
    assays = list(counts = Matrix::Matrix(0, nrow = 1, ncol = length(bc_all), sparse = TRUE)),
    colData = DataFrame(barcode = bc_all, in_tissue = as.logical(pos$in_tissue),
                        array_row = pos$array_row, array_col = pos$array_col, sum_umi = 0),
    spatialCoords = as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]))
  colnames(spe8) <- bc_all
  spe8$sample_id <- "s"

  metrics <- function(m) {
    s <- Matrix::colSums(m)
    data.frame(sum = s, det = Matrix::colSums(m > 0), pct = ifelse(s > 0, 100 * Matrix::colSums(m[mt, , drop = FALSE]) / s, 0))
  }
  if (bin == 8) {
    base <- metrics(m_all); skel <- spe8
  } else {
    skel <- aggregate_full_array(spe8, cd0, 16)
    key8 <- sprintf("s_016um_%05d_%05d-1", as.integer(pos$array_row %/% 2), as.integer(pos$array_col %/% 2))
    g <- match(key8, skel$barcode)
    M <- Matrix::sparseMatrix(i = seq_along(g), j = g, x = 1, dims = c(length(g), ncol(skel)))
    base <- metrics(m_all %*% M)
  }
  rm(m_all); invisible(gc())
  filelock::unlock(slot)
  saveRDS(list(base = base, skel = skel), paste0(cache, ".tmp"))
  file.rename(paste0(cache, ".tmp"), cache)
}
filelock::unlock(cell_lock)
off <- !as.logical(skel$in_tissue)

for (cond in todo) {
  t0 <- Sys.time()
  cd <- file.path(cell, cond)
  if (verify) {
    full <- read_full_array(cd, 8)
    if (bin != 8) full <- aggregate_full_array(full, cd, 16)
    stopifnot(identical(colnames(full), colnames(skel)))
    lg <- log10(full$sum_umi + 1)
    th <- isOutlier(lg, subset = as.logical(full$in_tissue), batch = full$sample_id, type = "lower", nmads = 3)
    orig <- lg < attr(th, "thresholds")["lower", "s"]
    hit <- spatial_steps(skel, orig, paste0(bin, "um"))
    if (bin != 8) hit <- project(hit, cd, bin)
    arm <- if (bin == 8) "19_spatialartifacts_8um" else "20_spatialartifacts_16um"
    ref <- readLines(file.path(cell, "results", arm, paste0(cond, ".txt")))
    say("VERIFY %s %s bin %d: copied steps reproduce arm %s: %s (%d against %d), %.0f s", basename(cell), cond, bin,
        substr(arm, 1, 2), setequal(hit, ref), length(hit), length(ref), as.numeric(difftime(Sys.time(), t0, units = "secs")))
    next
  }
  spe_in <- add_metrics(read_bin(cd, bin))
  nt <- !spe_in$is_truncated
  d <- perCellQCFilters(colData(spe_in)[nt, , drop = FALSE], sub.fields = "subsets_Mito_percent")
  thr <- c(lib = unname(attr(d$low_lib_size, "thresholds")["lower"]),
           det = unname(attr(d$low_n_features, "thresholds")["lower"]),
           mito = unname(attr(d$high_subsets_Mito_percent, "thresholds")["higher"]))
  met <- base
  i <- match(spe_in$barcode, skel$barcode)
  stopifnot(!anyNA(i), all(!off[i]))
  met$sum[i] <- spe_in$sum; met$det[i] <- spe_in$detected; met$pct[i] <- spe_in$subsets_Mito_percent
  flag <- unname(met$sum < thr["lib"] | met$det < thr["det"] | met$pct > thr["mito"])
  ours <- spe_in$barcode[nt][which(d$discard)]
  check1 <- setequal(ours, skel$barcode[flag & !off & skel$barcode %in% spe_in$barcode[nt]])
  hit <- spatial_steps(skel, flag, paste0(bin, "um"))
  if (bin != 8) hit <- project(hit, cd, bin)
  writeLines(hit, file.path(outdir, paste0(cond, ".txt")))
  say("%s %s bin %d: check1 %s, step 1 %d in tissue, final %d (8 um), %.0f s",
      basename(cell), cond, bin, check1, sum(flag & !off), length(hit),
      as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
