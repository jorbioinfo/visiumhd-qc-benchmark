suppressPackageStartupMessages(library(scuttle))

x  <- c(rep(0.1, 50), rep(5, 50))
it <- c(rep(0L, 50), rep(1L, 50))

th <- function(sub) attr(isOutlier(x, subset = sub, type = "lower", nmads = 3),
                         "thresholds")["lower"]

cat(sprintf("subset = integer 0/1        threshold %+.3f\n", th(it)))
cat(sprintf("subset = rep(1L, 100)       threshold %+.3f   <- identical\n", th(rep(1L, 100))))
cat(sprintf("subset = as.logical(0/1)    threshold %+.3f   <- intended\n", th(as.logical(it))))
cat(sprintf("\nx[1] = %.1f; the first two collapse onto it, so MAD = 0.\n", x[1]))

stopifnot(identical(unname(th(it)), unname(th(rep(1L, 100)))))
cat("\nFix: subset = as.logical(colData(spe)$in_tissue)\n")
