read_bin <- function(cond_dir, bin_um, keep_off_tissue = FALSE) {
  d <- file.path(cond_dir, sprintf("square_%03dum", bin_um))
  m <- Matrix::readMM(gzfile(file.path(d, "filtered_feature_bc_matrix", "matrix.mtx.gz")))
  bc <- readLines(gzfile(file.path(d, "filtered_feature_bc_matrix", "barcodes.tsv.gz")))
  ft <- read.delim(gzfile(file.path(d, "filtered_feature_bc_matrix", "features.tsv.gz")),
                   header = FALSE, col.names = c("id", "symbol", "type"))
  pos <- as.data.frame(arrow::read_parquet(file.path(d, "spatial", "tissue_positions.parquet")))
  pos <- pos[match(bc, pos$barcode), ]
  in_tis <- if ("in_tissue" %in% colnames(pos)) pos$in_tissue else rep(1L, nrow(pos))
  spe <- SpatialExperiment(
    assays = list(counts = as(m, "CsparseMatrix")),
    colData = DataFrame(barcode = bc, array_row = pos$array_row,
                        array_col = pos$array_col, is_truncated = pos$is_truncated,
                        in_tissue = in_tis,
                        n_children = if ("n_children" %in% colnames(pos))
                                       pos$n_children else rep(1L, nrow(pos))),
    rowData = DataFrame(ft),
    spatialCoords = as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]))
  colnames(spe) <- bc
  spe$sample_id <- "s"
  if (!keep_off_tissue) spe <- spe[, spe$in_tissue == 1]
  spe
}

add_metrics <- function(spe) {
  mito <- grepl("^MT-", rowData(spe)$symbol, ignore.case = TRUE)
  colData(spe) <- cbind(colData(spe), perCellQCMetrics(spe, subsets = list(Mito = mito)))
  spe$log_sum <- log1p(spe$sum)
  spe$subsets_Mito_percent[is.na(spe$subsets_Mito_percent)] <- 0
  spe$sum_umi  <- spe$sum
  spe$sum_gene <- spe$detected
  spe
}

drop_truncated <- function(spe) spe[, !spe$is_truncated]

run_mad <- function(spe) {
  spe <- add_metrics(drop_truncated(spe))
  d <- perCellQCFilters(colData(spe), sub.fields = "subsets_Mito_percent")
  spe$barcode[which(d$discard)]
}

run_spotsweeper <- function(spe, k = 36) {
  spe <- add_metrics(drop_truncated(spe))
  specs <- list(list(m = "sum",                  dir = "lower",  lg = TRUE),
                list(m = "detected",             dir = "lower",  lg = TRUE),
                list(m = "subsets_Mito_percent", dir = "higher", lg = FALSE))
  for (x in specs) {
    spe <- localOutliers(spe, metric = x$m, direction = x$dir, log = x$lg,
                         n_neighbors = k)
  }
  flags <- paste0(c("sum", "detected", "subsets_Mito_percent"), "_outliers")
  spe$discard <- Reduce(`|`, lapply(flags, function(f) as.logical(colData(spe)[[f]])))
  spe$barcode[which(spe$discard)]
}

check_source <- function(cond_dir) {
  a <- as.data.frame(arrow::read_parquet(
    file.path(SOURCE_8UM, "spatial", "tissue_positions.parquet")))
  a <- sort(a$barcode[a$in_tissue == 1])
  b <- sort(readLines(gzfile(file.path(cond_dir, "square_008um",
                                       "filtered_feature_bc_matrix", "barcodes.tsv.gz"))))
  if (!identical(a, b))
    stop(sprintf(paste("SOURCE_8UM does not match this benchmark: %d in-tissue bins",
                       "in the source vs %d in the condition (%s).",
                       "Set SOURCE_8UM to this sample's square_008um."),
                 length(a), length(b), SOURCE_8UM))
  invisible(TRUE)
}

read_full_array <- function(cond_dir, bin_um) {
  stopifnot(bin_um == 8)
  check_source(cond_dir)
  d <- file.path(cond_dir, "square_008um", "filtered_feature_bc_matrix")
  m_dmg <- Matrix::readMM(gzfile(file.path(d, "matrix.mtx.gz")))
  bc_dmg <- readLines(gzfile(file.path(d, "barcodes.tsv.gz")))
  umi_dmg <- Matrix::colSums(m_dmg); rm(m_dmg); invisible(gc())

  raw <- file.path(SOURCE_8UM, "raw_feature_bc_matrix")
  if (!file.exists(file.path(raw, "matrix.mtx.gz")))
    stop(sprintf(paste("this arm needs the whole capture area, but %s has no",
                       "raw_feature_bc_matrix. Run scripts/64_fetch_raw_8um.sh for",
                       "this dataset."), SOURCE_8UM))
  m_all <- Matrix::readMM(gzfile(file.path(raw, "matrix.mtx.gz")))
  bc_all <- readLines(gzfile(file.path(raw, "barcodes.tsv.gz")))
  umi <- Matrix::colSums(m_all); rm(m_all); invisible(gc())
  umi[match(bc_dmg, bc_all)] <- umi_dmg

  pos <- as.data.frame(arrow::read_parquet(
    file.path(SOURCE_8UM, "spatial", "tissue_positions.parquet")))
  pos <- pos[match(bc_all, pos$barcode), ]

  spe <- SpatialExperiment(
    assays = list(counts = Matrix::Matrix(0, nrow = 1, ncol = length(bc_all), sparse = TRUE)),
    colData = DataFrame(barcode = bc_all, in_tissue = as.logical(pos$in_tissue),
                        array_row = pos$array_row, array_col = pos$array_col,
                        sum_umi = umi),
    spatialCoords = as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]))
  colnames(spe) <- bc_all
  spe$sample_id <- "s"
  spe
}

aggregate_full_array <- function(spe8, cond_dir, bin_um) {
  f <- bin_um %/% 8
  r <- spe8$array_row %/% f; cc <- spe8$array_col %/% f
  key <- r * 100000 + cc
  u <- unique(key); g <- match(key, u)
  ur <- as.integer(u %/% 100000); uc <- as.integer(u %% 100000)
  n <- tabulate(g, length(u))
  tis <- as.vector(rowsum(as.integer(spe8$in_tissue), g)) > 0
  bc <- sprintf("s_%03dum_%05d_%05d-1", as.integer(bin_um), ur, uc)
  valid <- readLines(gzfile(file.path(cond_dir, sprintf("square_%03dum", bin_um),
                                      "filtered_feature_bc_matrix", "barcodes.tsv.gz")))
  if (!setequal(bc[tis], valid))
    stop(sprintf("%d um tissue bins from the full array (%d) differ from %s (%d)",
                 bin_um, sum(tis), sprintf("square_%03dum", bin_um), length(valid)))
  xy <- rowsum(spatialCoords(spe8), g) / n
  rownames(xy) <- NULL
  spe <- SpatialExperiment(
    assays = list(counts = Matrix::Matrix(0, nrow = 1, ncol = length(u), sparse = TRUE)),
    colData = DataFrame(barcode = bc, in_tissue = tis, array_row = ur, array_col = uc,
                        sum_umi = as.vector(rowsum(spe8$sum_umi, g))),
    spatialCoords = xy)
  colnames(spe) <- bc
  spe$sample_id <- "s"
  spe
}

run_spatialartifacts_v <- function(spe, bin_um) {
  stopifnot(bin_um %in% c(8, 16))
  stopifnot(all(startsWith(spe$barcode, sprintf("s_%03dum_", bin_um))))
  det <- detectEdgeArtifacts(spe, platform = "visiumhd",
                             resolution = paste0(bin_um, "um"),
                             qc_metric = "sum_umi", samples = "sample_id",
                             batch_var = "sample_id", mad_threshold = 3,
                             name = "edge_artifact", verbose = FALSE)
  cls <- classifyEdgeArtifacts(det, min_spots = 20, name = "edge_artifact")
  v <- colData(cls)[["edge_artifact_classification"]]
  flag <- !is.na(v) & v != "not_artifact"
  cls$barcode[which(flag & as.logical(cls$in_tissue))]
}

project <- function(bad_parents, cond_dir, bin_um) {
  if (length(bad_parents) == 0) return(character(0))
  stopifnot(all(startsWith(bad_parents, sprintf("s_%03dum_", bin_um))))
  d_hi <- file.path(cond_dir, "square_008um")
  valid <- readLines(gzfile(file.path(d_hi, "filtered_feature_bc_matrix", "barcodes.tsv.gz")))
  barcode_based_deconvolve(bad_parents, low_res = bin_um, high_res = 8,
                           valid_barcodes = valid)
}

barcode_based_deconvolve <- function(barcode_list, low_res, high_res,
                                     valid_barcodes = NULL) {

  if (missing(barcode_list) || !is.vector(barcode_list) || !is.character(barcode_list)) {
    stop("`barcode_list` must be a non-empty character vector.")
  }
  if (length(barcode_list) == 0) {
    return(character(0))
  }
  if (!is.numeric(low_res) || !is.numeric(high_res) ||
      length(low_res) != 1 || length(high_res) != 1 ||
      low_res != floor(low_res) || high_res != floor(high_res) ||
      low_res <= 0 || high_res <= 0) {
    stop("`low_res` and `high_res` must be positive integer values.")
  }

  if (low_res < high_res) {
    stop("low_res must be larger than high_res.")
  }
  if (low_res %% high_res != 0) {
    stop("low_res must be a multiple of high_res.")
  }
  factor <- low_res / high_res
  prefix <- sprintf("%03dum", high_res)

  if (any(is.na(barcode_list)) || !all(grepl("_", barcode_list, fixed = TRUE))) {
    stop("Each element of `barcode_list` must be a non-NA string containing '_' separators.")
  }

  parts <- strsplit(sub("-1$", "", barcode_list), "_", fixed = TRUE)
  if (any(lengths(parts) < 4)) {
    stop(paste0("Barcode '", barcode_list[which(lengths(parts) < 4)[1]],
                "' does not conform to expected pattern."))
  }
  row_low <- suppressWarnings(as.integer(vapply(parts, `[`, character(1), 3)))
  col_low <- suppressWarnings(as.integer(vapply(parts, `[`, character(1), 4)))
  if (anyNA(row_low) || anyNA(col_low)) {
    stop(paste0("Failed to parse numeric row/col from barcode '",
                barcode_list[which(is.na(row_low) | is.na(col_low))[1]], "'."))
  }

  offs <- seq_len(factor) - 1L
  grid <- expand.grid(off_c = offs, off_r = offs, i = seq_along(row_low))
  discard_barcodes <- sprintf("s_%s_%05d_%05d-1", prefix,
                              row_low[grid$i] * factor + grid$off_r,
                              col_low[grid$i] * factor + grid$off_c)
  discard_barcodes <- unique(discard_barcodes)

  if (!is.null(valid_barcodes)) {
    n_before <- length(discard_barcodes)
    discard_barcodes <- discard_barcodes[discard_barcodes %in% valid_barcodes]
    n_dropped <- n_before - length(discard_barcodes)
    if (n_dropped > 0) {
      message(sprintf("barcode_based_deconvolve: dropped %d generated barcode(s) absent from `valid_barcodes`.",
                      n_dropped))
    }
  }
  return(discard_barcodes)
}

