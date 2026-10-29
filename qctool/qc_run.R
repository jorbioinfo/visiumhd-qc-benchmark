suppressPackageStartupMessages({
  library(Matrix); library(arrow); library(jsonlite)
  library(SpotSweeper); library(SpatialExperiment); library(scuttle); library(SpatialArtifacts)
})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d = NULL) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "qc_core.R"))

VIEWS <- get_arg("--views"); TEST <- get_arg("--test"); OUT <- get_arg("--out")
BIN <- as.integer(get_arg("--bin", "8")); K <- as.integer(get_arg("--k", "48"))
SOURCE_8UM <- get_arg("--source8", file.path(VIEWS, "square_008um"))
stopifnot(!is.null(VIEWS), !is.null(TEST), !is.null(OUT))

removed <- switch(TEST,
  mad = if (BIN == 8) run_mad(read_bin(VIEWS, 8))
        else project(run_mad(read_bin(VIEWS, BIN)), VIEWS, BIN),
  spotsweeper = if (BIN == 8) run_spotsweeper(read_bin(VIEWS, 8), k = K)
                else project(run_spotsweeper(read_bin(VIEWS, BIN), k = K), VIEWS, BIN),
  spatialartifacts = if (BIN == 8) run_spatialartifacts_v(read_full_array(VIEWS, 8), 8)
                     else if (BIN == 16) project(run_spatialartifacts_v(
                       aggregate_full_array(read_full_array(VIEWS, 8), VIEWS, 16), 16), VIEWS, 16)
                     else stop("SpatialArtifacts accepts 8 or 16 um bins only"),
  stop("unknown --test ", TEST))
writeLines(sort(unique(removed)), OUT)
cat(sprintf("%s %d um%s: %d bins removed\n", TEST, BIN,
            if (TEST == "spotsweeper") sprintf(", %d neighbours", K) else "", length(unique(removed))))
