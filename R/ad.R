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
