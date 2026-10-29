suppressPackageStartupMessages({
  library(SingleCellExperiment); library(miQC); library(flexmix)
})
set.seed(1)
n <- 2000
cd <- DataFrame(
  detected = rpois(n, 400),
  subsets_Mito_percent = c(rnorm(n * 0.8, 5, 2), rnorm(n * 0.2, 40, 8))
)
sce <- SingleCellExperiment(colData = cd)

cat("--- as documented: pass the column name explicitly ---\n")
print(try(mixtureModel(sce, subsets_mito_percent = "subsets_Mito_percent"), silent = TRUE))

cat("\n--- what actually works: rename the column ---\n")
sce$subsets_mito_percent <- sce$subsets_Mito_percent
m <- mixtureModel(sce)
cat("components:", length(m@components), "\n")
cat("parameters:\n"); print(parameters(m))
