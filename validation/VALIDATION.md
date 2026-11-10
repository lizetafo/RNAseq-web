# Validation

**Question:** does the app give exactly the same results as running DESeq2 directly?

For each example dataset, [`compare_to_deseq2.R`](compare_to_deseq2.R) runs:
- **(a)** the app's own functions (`read_file`, `qc_check`, `prepare_meta`, `run_deseq`, `extract_results`), loaded from `app.R`, and
- **(b)** a plain DESeq2 analysis written independently: `DESeqDataSetFromMatrix` → filter genes with < 10 total counts → `DESeq` → `results(contrast, alpha = 0.05)`.

It then compares every result column (`baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`) and the normalized counts.

Re-run it in RStudio from the repository folder:
```r
source("validation/compare_to_deseq2.R")
```

## Results (R 4.3, DESeq2 1.42)

| Dataset | Comparison | Genes tested | Same genes | Largest difference: results / normalized counts | Significant genes (padj < 0.05): app / DESeq2 |
|---|---|---|---|---|---|
| prostate_hypoxia (`~Condition`) | Hypoxia vs Normoxia | 25,062 | ✅ | 0 / 0 | 302 / 302 |
| prostate_hypoxia (`~CellLine + Condition`) | Hypoxia vs Normoxia | 25,062 | ✅ | 0 / 0 | 4,895 / 4,895 |
| tumor_vs_normal | cancer vs normal | 28,225 | ✅ | 0 / 0 | 6,834 / 6,834 |

The app's results are **identical** to plain DESeq2 for both datasets and for a two-factor design. The app adds checks, annotation and plots without changing the statistics.

## Full app run
The app was also run in a browser: upload → choose control group → Run → all plots → download of the tables.
- After upload, the **control group is left blank**. Clicking Run before choosing one shows *"Please choose the reference (control) group"*.
- Choosing **Normoxia** filled in the contrast `Condition, Hypoxia, Normoxia`.
- The downloaded results matched plain DESeq2 to within 10⁻¹⁴ (rounding from saving the CSV), with the same 302 significant genes.
- A sample sheet with a stray space in a group label (`Activated ` vs `Activated`) was read as **two** groups, not three, because labels are trimmed on upload.
- Gene symbols were assigned to about 76% of Ensembl IDs. IDs without a symbol are mostly non-coding or novel genes.
- **Quality-control and extra plots:** every new plot (PCA, sample distances, library size, dispersion, p-value histogram, MA plot, Gene Explorer) and its download were checked on both datasets, including the 37-sample colorectal dataset (analysis ≈ 40 s, rlog heatmap ≈ 1 min extra). The **Load example data** button loads the prostate dataset.
