# Methods

## Where the app fits in an RNA-seq workflow

```
FASTQ reads
   │  FastQC         quality check of raw reads
   │  Trimmomatic    trimming of adapters and low-quality bases
   │  HISAT2         alignment to the reference genome (SAM/BAM)
   │  featureCounts  gene-level read counts (Subread)
   ▼
Count matrix  ──►  THIS APP: QC checks → DESeq2 → annotation → plots and tables
```

Read processing up to the count matrix is done on the command line (Linux). The app starts from the **raw gene count matrix** plus a sample sheet.

## Analysis steps

0. **Getting the data.** Files are read with `read_counts`, which handles single tables, featureCounts output, per-sample HTSeq/STAR files and archives, and refuses non-integer (normalized) or negative values (`check_count_values`). Alternatively, `load_geo` downloads NCBI's raw count matrix for a GEO series (human GRCh38.p13 or mouse GRCm39, Entrez gene IDs) and builds the sample sheet from the Series Matrix sample characteristics (`parse_series_matrix`).
1. **Reading the data** (`read_file`). Counts and metadata are read from CSV, TSV or TXT. Gene IDs become row names, and counts are stored as integers. Missing values are set to zero if requested; otherwise the app stops. Duplicated sample names stop the analysis, and duplicated gene IDs are removed with a message.
2. **Quality checks** (`qc_check`). The count columns and metadata rows must contain the same samples, in the same order.
3. **Preparing the design** (`prepare_meta`). The grouping columns are converted to factors, and the reference (control) level is set.
4. **Differential expression** (`run_deseq`). A `DESeqDataSet` is built with the design formula. Genes with fewer than 10 counts in total across all samples are removed, then `DESeq()` estimates size factors and dispersions and fits a negative binomial model.
5. **Results** (`extract_results`). `results()` is called for the chosen contrast with the chosen alpha, and normalized counts are extracted. Adjusted p-values use Benjamini–Hochberg.
6. **Annotation** (`annotate_results`). Gene IDs are mapped to gene symbols with `AnnotationDbi::mapIds()` and the selected organism database. If the ID type is not given, it is guessed from the ID format.
7. **Top genes** (`top_genes`). Genes are sorted by adjusted p-value.
8. **Quality control.** Data are variance-stabilized with `vst()` (blind to the design), falling back to `varianceStabilizingTransformation()` for very small datasets. PCA uses `prcomp` on the most variable genes, as `DESeq2::plotPCA` does. Sample distances are Euclidean distances on VST values. Dispersion estimates come from `plotDispEsts()`.
9. **Optional settings.** A subset of samples can be chosen before the analysis (`subset_samples`); design factors left with one value are removed (`drop_single_level_terms`). Groups with fewer than 3 samples trigger a warning (`replicate_warning`). Fold changes can be shrunk with `lfcShrink(type = "normal")` (`shrink_results`), which leaves p-values unchanged. Any DESeq2 coefficient can be shown with `results(name = ...)`.
10. **Pathway analysis** (`run_gsea`). Genes are ranked by the DESeq2 Wald statistic. Where several IDs share a gene symbol, the one with the largest absolute statistic is kept. `fgsea()` tests the 50 MSigDB Hallmark gene sets (size 15–500).
11. **Plots.** The volcano and MA plots use ggplot2 and ggrepel. The heatmap uses ComplexHeatmap on z-scored normalized, VST or rlog values; rlog is computed only when selected, because it is slow. The bar plot shows log2 fold changes. The Gene Explorer plots normalized counts for one gene.

## Software
R · Shiny · DESeq2 · fgsea · AnnotationDbi · org.\*.eg.db · ComplexHeatmap · circlize · ggplot2 · ggrepel · dplyr · DT

