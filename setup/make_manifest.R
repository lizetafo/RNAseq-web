# ================================================================
# Creates manifest.json for hosting the app on Posit Connect Cloud.
#
# Run ONCE on your own computer, in RStudio, from the repository folder:
#   source("setup/make_manifest.R")
# Then upload the new manifest.json to GitHub (same folder as app.R).
# Re-run it only if you add new packages to the app.
# ================================================================

if (!requireNamespace("rsconnect", quietly = TRUE)) install.packages("rsconnect")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Every package the app can use must be installed here so it is recorded
cran <- c("shiny", "dplyr", "ggplot2", "ggrepel", "forcats", "RColorBrewer", "circlize", "DT", "readxl")
bioc <- c("DESeq2", "AnnotationDbi", "ComplexHeatmap", "fgsea",
          "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
          "org.Dm.eg.db", "org.Dr.eg.db", "org.Sc.sgd.db")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

# Tell rsconnect where Bioconductor packages come from
options(repos = BiocManager::repositories())

rsconnect::writeManifest(
  appDir = ".",
  appFiles = c("app.R",
               "example_data/prostate_hypoxia/counts.csv",
               "example_data/prostate_hypoxia/metadata.csv",
               "gene_sets/hallmark_human.csv")
)
cat("\nmanifest.json created. Upload it to GitHub, in the same folder as app.R.\n")
