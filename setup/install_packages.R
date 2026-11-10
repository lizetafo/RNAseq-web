# ================================================================
# One-time installation of the packages the app uses.
# Run once in the RStudio Console, from the repository folder:
#   source("setup/install_packages.R")
# ================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c("shiny", "dplyr", "ggplot2", "ggrepel", "forcats",
                   "RColorBrewer", "circlize", "DT"))

# DESeq2, annotation and heatmap packages, plus the human gene database.
BiocManager::install(c("DESeq2", "AnnotationDbi", "ComplexHeatmap", "org.Hs.eg.db"),
                     ask = FALSE, update = FALSE)

# Other organisms are installed automatically the first time you select
# them in the app. To install them in advance, uncomment what you need:
# BiocManager::install("org.Mm.eg.db")   # mouse
# BiocManager::install("org.Rn.eg.db")   # rat
# BiocManager::install("org.Dm.eg.db")   # fruit fly
# BiocManager::install("org.Dr.eg.db")   # zebrafish
# BiocManager::install("org.Sc.sgd.db")  # yeast

cat("\nDone. Restart R (Session > Restart R), open app.R and click Run App.\n")
