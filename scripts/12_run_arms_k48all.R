suppressPackageStartupMessages({
  library(Matrix); library(arrow); library(jsonlite); library(parallel)
  library(SpotSweeper); library(SpatialExperiment); library(scuttle)
  library(SpatialArtifacts)
  library(miQC); library(flexmix)
})
source(file.path(Sys.getenv("VHDQC_ROOT", "/work"), "R", "ad.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }
BENCH   <- get_arg("--benchmark", "/work/benchmark")
RESULTS <- get_arg("--results",   file.path(BENCH, "results"))
WORKERS <- as.integer(get_arg("--workers", "8"))
ONLY    <- get_arg("--arms", "")
CONDS   <- get_arg("--conds", "")
HAVE_SCRAPPER <- requireNamespace("scrapper", quietly = TRUE)
LOGDIR  <- file.path(RESULTS, "logs")
dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)

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

run_findartifacts <- function(spe) {
  spe <- add_metrics(drop_truncated(spe))
  spe$expr_chrM       <- spe$subsets_Mito_sum
  spe$expr_chrM_ratio <- spe$subsets_Mito_percent / 100
  spe <- tryCatch(
    findArtifacts(spe, mito_percent = "expr_chrM_ratio", mito_sum = "expr_chrM",
                  samples = "sample_id", n_order = 5, shape = "square",
                  log = TRUE, name = "artifact"),
    error = function(e) {
      message("findArtifacts: log=TRUE failed (", conditionMessage(e),
              "); retrying with log=FALSE")
      findArtifacts(spe, mito_percent = "expr_chrM_ratio", mito_sum = "expr_chrM",
                    samples = "sample_id", n_order = 5, shape = "square",
                    log = FALSE, name = "artifact")
    })
  v <- colData(spe)[["artifact"]]
  flag <- if (is.logical(v)) v else as.character(v) %in% c("TRUE", "artifact", "1")
  spe$barcode[which(flag)]
}

run_scrapper <- function(spe) {
  if (!HAVE_SCRAPPER) stop("scrapper not installed in this image")
  spe <- drop_truncated(spe)
  mito <- grepl("^MT-", rowData(spe)$symbol, ignore.case = TRUE)
  m <- scrapper::computeRnaQcMetrics(counts(spe), subsets = list(Mito = which(mito)))
  th <- scrapper::suggestRnaQcThresholds(m)
  keep <- scrapper::filterRnaQcMetrics(th, m)
  spe$barcode[!keep]
}

SOURCE_8UM <- local({
  env <- Sys.getenv("SOURCE_8UM", "")
  if (nzchar(env)) return(env)
  f <- file.path(BENCH, "injection_params.json")
  if (file.exists(f)) {
    src <- jsonlite::fromJSON(f)$src
    if (!is.null(src) && dir.exists(file.path(src, "square_008um")))
      return(file.path(src, "square_008um"))
  }
  "/work/runs/armA_sameinput_bin64/outs/binned_outputs/square_008um"
})

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

drop_truncated <- function(spe) spe[, !spe$is_truncated]

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

run_miqc <- function(spe) {
  spe <- add_metrics(drop_truncated(spe))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = SummarizedExperiment::assay(spe, "counts")),
    colData = colData(spe))
  colnames(sce) <- spe$barcode
  sce$subsets_mito_percent <- sce$subsets_Mito_percent

  model <- withCallingHandlers(
    miQC::mixtureModel(sce, model_type = "linear"),
    warning = function(w) {
      message("  miQC warning: ", conditionMessage(w)); invokeRestart("muffleWarning")
    })

  if (length(model@components) < 2L) {
    message("  miQC: flexmix returned ", length(model@components),
            " component(s); two-distribution assumption not met, nothing filtered")
    return(character(0))
  }

  kept <- tryCatch(
    miQC::filterCells(sce, model, posterior_cutoff = 0.75,
                      keep_all_below_boundary = TRUE, enforce_left_cutoff = TRUE,
                      verbose = FALSE),
    error = function(e) {
      message("  miQC: filterCells declined to filter (", conditionMessage(e), ")")
      NULL
    })
  if (is.null(kept)) return(character(0))
  setdiff(colnames(sce), colnames(kept))
}

run_fixed_umi <- function(spe, bin_um = 8) {
  spe <- add_metrics(drop_truncated(spe))
  spe$barcode[spe$sum <= 100 * (bin_um / 8)^2]
}

run_mad_percovered <- function(spe) {
  spe <- add_metrics(drop_truncated(spe))
  nc <- pmax(as.numeric(spe$n_children), 1)
  cd <- colData(spe)
  cd$sum <- cd$sum / nc
  cd$detected <- cd$detected / nc
  d <- perCellQCFilters(cd, sub.fields = "subsets_Mito_percent")
  spe$barcode[which(d$discard)]
}

run_mad_excludepartial <- function(spe, full_children) {
  spe <- add_metrics(drop_truncated(spe))
  full <- as.numeric(spe$n_children) >= full_children
  if (!any(full)) return(character(0))
  d <- perCellQCFilters(colData(spe)[full, , drop = FALSE],
                        sub.fields = "subsets_Mito_percent")
  spe$barcode[which(full)][which(d$discard)]
}

run_mad <- function(spe) {
  spe <- add_metrics(drop_truncated(spe))
  d <- perCellQCFilters(colData(spe), sub.fields = "subsets_Mito_percent")
  spe$barcode[which(d$discard)]
}

run_spatialartifacts <- function(spe, bin_um, fix_in_tissue = FALSE) {
  stopifnot(bin_um %in% c(8, 16))
  spe <- add_metrics(drop_truncated(spe))
  if (fix_in_tissue) spe$in_tissue <- as.logical(spe$in_tissue)
  spe <- detectEdgeArtifacts(spe, platform = "visiumhd",
                             resolution = paste0(bin_um, "um"),
                             qc_metric = "sum_gene", verbose = FALSE)
  flag <- as.logical(colData(spe)[["edge_artifact"]])
  spe$barcode[which(flag & spe$in_tissue == 1)]
}

project <- function(bad_parents, cond_dir, bin_um) {
  if (length(bad_parents) == 0) return(character(0))
  stopifnot(all(startsWith(bad_parents, sprintf("s_%03dum_", bin_um))))
  d_hi <- file.path(cond_dir, "square_008um")
  valid <- readLines(gzfile(file.path(d_hi, "filtered_feature_bc_matrix", "barcodes.tsv.gz")))
  barcode_based_deconvolve(bad_parents, low_res = bin_um, high_res = 8,
                           valid_barcodes = valid)
}

arms <- list(
  "01_no_qc"            = function(cd) character(0),
  "02_mad_8um"          = function(cd) run_mad(read_bin(cd, 8)),
  "03_spotsweeper_8um"  = function(cd) run_spotsweeper(read_bin(cd, 8)),
  "04_spotsweeper_16um" = function(cd) project(run_spotsweeper(read_bin(cd, 16)), cd, 16),
  "05_spotsweeper_32um" = function(cd) project(run_spotsweeper(read_bin(cd, 32)), cd, 32),
  "06_spotsweeper_64um" = function(cd) project(run_spotsweeper(read_bin(cd, 64)), cd, 64),
  "07_spatialartifacts" = function(cd) run_spatialartifacts(read_bin(cd, 8, keep_off_tissue = TRUE), 8),

  "24_spotsweeper_8um_k48"   = function(cd) run_spotsweeper(read_bin(cd, 8), k = 48),
  "21_spotsweeper_8um_k144"  = function(cd) run_spotsweeper(read_bin(cd, 8), k = 144),
  "22_spotsweeper_8um_k576"  = function(cd) run_spotsweeper(read_bin(cd, 8), k = 576),
  "23_spotsweeper_8um_k2304" = function(cd) run_spotsweeper(read_bin(cd, 8), k = 2304),

  "25_spotsweeper_16um_k48"  = function(cd) project(run_spotsweeper(read_bin(cd, 16), k = 48), cd, 16),
  "26_spotsweeper_32um_k48"  = function(cd) project(run_spotsweeper(read_bin(cd, 32), k = 48), cd, 32),
  "27_spotsweeper_64um_k48"  = function(cd) project(run_spotsweeper(read_bin(cd, 64), k = 48), cd, 64),
  "28_spotsweeper_128um_k48" = function(cd) project(run_spotsweeper(read_bin(cd, 128), k = 48), cd, 128),
  "44_spotsweeper_8um_k192"  = function(cd) run_spotsweeper(read_bin(cd, 8), k = 192),
  "45_spotsweeper_8um_k768"  = function(cd) run_spotsweeper(read_bin(cd, 8), k = 768),
  "46_spotsweeper_8um_k3072" = function(cd) run_spotsweeper(read_bin(cd, 8), k = 3072),

  "30_miqc_8um"         = function(cd) run_miqc(read_bin(cd, 8)),
  "31_miqc_64um"        = function(cd) project(run_miqc(read_bin(cd, 64)), cd, 64),
  "32_fixedumi_8um"     = function(cd) run_fixed_umi(read_bin(cd, 8), 8),
  "40_spotsweeper_128um" = function(cd) project(run_spotsweeper(read_bin(cd, 128)), cd, 128),
  "41_mad_128um"         = function(cd) project(run_mad(read_bin(cd, 128)), cd, 128),
  "34_mad_64um_percovered"    = function(cd) project(run_mad_percovered(read_bin(cd, 64)), cd, 64),
  "35_mad_64um_excludepartial" = function(cd) project(
      run_mad_excludepartial(read_bin(cd, 64), 64), cd, 64),
  "33_fixedumi_64um"    = function(cd) project(run_fixed_umi(read_bin(cd, 64), 64), cd, 64),

  "08_mad_16um"         = function(cd) project(run_mad(read_bin(cd, 16)), cd, 16),
  "09_mad_32um"         = function(cd) project(run_mad(read_bin(cd, 32)), cd, 32),
  "10_mad_64um"         = function(cd) project(run_mad(read_bin(cd, 64)), cd, 64),

  "16_spatialartifacts_8um_fixed" = function(cd)
      run_spatialartifacts(read_bin(cd, 8, keep_off_tissue = TRUE), 8, fix_in_tissue = TRUE),
  "11_findartifacts_8um"  = function(cd) run_findartifacts(read_bin(cd, 8)),
  "12_findartifacts_64um" = function(cd) project(run_findartifacts(read_bin(cd, 64)), cd, 64),
  "13_scrapper_8um"       = function(cd) run_scrapper(read_bin(cd, 8)),
  "19_spatialartifacts_8um"  = function(cd) run_spatialartifacts_v(read_full_array(cd, 8), 8),
  "20_spatialartifacts_16um" = function(cd) project(
      run_spatialartifacts_v(aggregate_full_array(read_full_array(cd, 8), cd, 16), 16), cd, 16),
  "17_scrapper_16um"      = function(cd) project(run_scrapper(read_bin(cd, 16)), cd, 16),
  "18_scrapper_32um"      = function(cd) project(run_scrapper(read_bin(cd, 32)), cd, 32),
  "14_scrapper_64um"      = function(cd) project(run_scrapper(read_bin(cd, 64)), cd, 64),
  "15_spatialartifacts_16um" = function(cd)
      project(run_spatialartifacts(read_bin(cd, 16, keep_off_tissue = TRUE), 16), cd, 16)
)

conds <- basename(list.dirs(BENCH, recursive = FALSE))
conds <- setdiff(conds, c("results", "scored"))
if (nzchar(CONDS)) {
  keepc <- trimws(strsplit(CONDS, ",")[[1]])
  missing_c <- setdiff(keepc, conds)
  if (length(missing_c)) stop("conditions not in this cell: ", paste(missing_c, collapse = ", "))
  conds <- keepc
  message(sprintf("conditions restricted to: %s", paste(conds, collapse = ", ")))
}

if (nzchar(ONLY)) {
  keep <- trimws(strsplit(ONLY, ",")[[1]])
  stopifnot(all(keep %in% names(arms)))
  arms <- arms[keep]
  message(sprintf("restricted to: %s", paste(keep, collapse = ", ")))
}

tasks <- expand.grid(arm = names(arms), cond = conds, stringsAsFactors = FALSE)
tasks$out <- file.path(RESULTS, tasks$arm, paste0(tasks$cond, ".txt"))
for (a in names(arms)) dir.create(file.path(RESULTS, a), showWarnings = FALSE, recursive = TRUE)
todo <- tasks[!file.exists(tasks$out), ]

cost <- c("23_spotsweeper_8um_k2304" = 6, "22_spotsweeper_8um_k576" = 5,
          "21_spotsweeper_8um_k144" = 4, "24_spotsweeper_8um_k48" = 4,
          "03_spotsweeper_8um" = 4, "07_spatialartifacts" = 3, "02_mad_8um" = 2,
          "08_mad_16um" = 2, "09_mad_32um" = 1, "10_mad_64um" = 1,
          "11_findartifacts_8um" = 4, "12_findartifacts_64um" = 1,
          "13_scrapper_8um" = 2, "14_scrapper_64um" = 1,
          "17_scrapper_16um" = 2, "18_scrapper_32um" = 1,
          "19_spatialartifacts_8um" = 0, "20_spatialartifacts_16um" = 0,
          "15_spatialartifacts_16um" = 2, "16_spatialartifacts_8um_fixed" = 3,
          "04_spotsweeper_16um" = 2, "05_spotsweeper_32um" = 1,
          "06_spotsweeper_64um" = 1, "01_no_qc" = 0,
          "30_miqc_8um" = 2, "31_miqc_64um" = 1,
          "32_fixedumi_8um" = 1, "33_fixedumi_64um" = 1,
          "34_mad_64um_percovered" = 1, "35_mad_64um_excludepartial" = 1,
          "40_spotsweeper_128um" = 1, "41_mad_128um" = 0,
          "46_spotsweeper_8um_k3072" = 7, "45_spotsweeper_8um_k768" = 5,
          "44_spotsweeper_8um_k192" = 4, "24_spotsweeper_8um_k48" = 4,
          "25_spotsweeper_16um_k48" = 2, "26_spotsweeper_32um_k48" = 1,
          "27_spotsweeper_64um_k48" = 1, "28_spotsweeper_128um_k48" = 0)
todo <- todo[order(-cost[todo$arm]), ]

message(sprintf("%d task(s) over %d worker(s); %d already done",
                nrow(todo), WORKERS, nrow(tasks) - nrow(todo)))

run_one <- function(i) {
  arm <- todo$arm[i]; cond <- todo$cond[i]; out <- todo$out[i]
  lf <- file.path(LOGDIR, sprintf("%s__%s.log", arm, cond))
  t0 <- Sys.time()
  res <- tryCatch({
    bad <- arms[[arm]](file.path(BENCH, cond))
    writeLines(bad, out)
    sprintf("%-22s %-16s %7d discarded  (%.1fs)", arm, cond, length(bad),
            as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }, error = function(e)
    sprintf("%-22s %-16s ERROR: %s", arm, cond, conditionMessage(e)))
  writeLines(res, lf)
  res
}

out <- mclapply(seq_len(nrow(todo)), run_one, mc.cores = WORKERS, mc.preschedule = FALSE)
for (l in out) message(l)
message("done")
