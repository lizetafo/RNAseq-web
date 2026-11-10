# ================================================================
# Validation: does the app give the same results as DESeq2 run directly?
#
# For each example dataset, this script runs
#   (a) the app's own functions (read_file, qc_check, prepare_meta,
#       run_deseq, extract_results), loaded from app.R, and
#   (b) a plain DESeq2 analysis written independently,
# then compares every result column and the normalized counts.
#
# Run in RStudio from the repository folder:
#   source("validation/compare_to_deseq2.R")
# ================================================================

suppressPackageStartupMessages(library(DESeq2))

# Load only the function definitions from app.R (the app is not launched)
app_env <- new.env()
for (expr in parse("app.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, app_env)
  }
}

datasets <- list(
  list(name = "prostate_hypoxia",     factor = "Condition",  case = "Hypoxia", control = "Normoxia"),
  list(name = "prostate_hypoxia",     factor = "Condition",  case = "Hypoxia", control = "Normoxia",
       design = "~CellLine + Condition", groups = c("Condition", "CellLine")),
  list(name = "tumor_vs_normal",      factor = "type",       case = "cancer",  control = "normal")
)

summary_rows <- list()
for (d in datasets) {
  counts_file <- file.path("example_data", d$name, "counts.csv")
  meta_file <- file.path("example_data", d$name, "metadata.csv")
  cat("\n==", d$name, "==\n")

  # (a) the app's functions
  dat <- app_env$read_file(counts_file, meta_file, na_to_zero = TRUE)
  app_env$qc_check(dat$counts, dat$metadata)
  ref_levels <- list(); ref_levels[[d$factor]] <- d$control
  design <- if (is.null(d$design)) paste0("~", d$factor) else d$design
  groups <- if (is.null(d$groups)) d$factor else d$groups
  meta <- app_env$prepare_meta(dat$metadata, Exp_groups = groups, ref_levels = ref_levels)
  dds_app <- suppressMessages(app_env$run_deseq(dat$counts, meta, Exp_design = as.formula(design)))
  app <- app_env$extract_results(dds_app, contrast = c(d$factor, d$case, d$control), alpha = 0.05)
  app_res <- as.data.frame(app$count_results)

  # (b) plain DESeq2, written independently
  cts <- as.matrix(read.csv(counts_file, check.names = FALSE, row.names = 1))
  md <- read.csv(meta_file, check.names = FALSE, row.names = 1)
  for (g in groups) md[[g]] <- factor(md[[g]])
  md[[d$factor]] <- relevel(md[[d$factor]], ref = d$control)
  dds <- DESeqDataSetFromMatrix(cts, md, design = as.formula(design))
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)
  ref_res <- as.data.frame(results(dds, contrast = c(d$factor, d$case, d$control), alpha = 0.05))

  same_genes <- identical(sort(rownames(ref_res)), sort(rownames(app_res)))
  g <- rownames(ref_res)
  max_diff <- max(sapply(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                         function(col) max(abs(ref_res[g, col] - app_res[g, col]), na.rm = TRUE)))
  norm_diff <- max(abs(counts(dds, normalized = TRUE)[g, ] - app$normalized_counts[g, colnames(dds)]))
  row <- data.frame(dataset = d$name, design = design, comparison = paste(d$case, "vs", d$control),
                    genes_tested = nrow(ref_res), same_genes = same_genes,
                    max_difference_results = signif(max_diff, 2),
                    max_difference_normalized_counts = signif(norm_diff, 2),
                    significant_app = sum(app_res$padj < 0.05, na.rm = TRUE),
                    significant_deseq2 = sum(ref_res$padj < 0.05, na.rm = TRUE))
  print(row, row.names = FALSE)
  summary_rows[[length(summary_rows) + 1]] <- row
}

summary <- do.call(rbind, summary_rows)
utils::write.csv(summary, file.path("validation", "validation_results.csv"), row.names = FALSE)
cat("\nIdentical results for all datasets:",
    all(summary$same_genes & summary$max_difference_results < 1e-8 &
          summary$significant_app == summary$significant_deseq2), "\n")
