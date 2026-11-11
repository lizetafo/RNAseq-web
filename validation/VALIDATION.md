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
- **Subsets, shrinkage, advanced results and pathway analysis:** checked in the browser. Keeping only PC3 samples with design `~CellLine + Condition` correctly removed `CellLine` (design used: `~Condition`), warned about 2 replicates per group, and gave the same 3,689 significant genes as DESeq2 run on the PC3 samples directly. Shrinkage changed fold changes but left p-values identical. With design `~CellLine + Condition + CellLine:Condition`, the interaction term `CellLinePC3.ConditionHypoxia` could be selected, and the app correctly reported that shrinkage is not available for interaction designs.
- **Pathway analysis as a biology check:** in PC3 hypoxia vs normoxia, the top Hallmark pathway is **HYPOXIA** (NES 2.72, padj 3.3e-45), followed by glycolysis and mTORC1 signaling, as expected for cells grown in low oxygen. Hypoxia genes named in the original study (ADM, ANGPTL4, CA9) are strongly up (log2 fold change 4.9–6.6).
- **Input formats:** featureCounts output, a GEO-style `_RAW.tar` of gzipped per-sample HTSeq files, and several separate per-sample files all gave count tables **identical** to the original matrix (HTSeq summary rows such as `__no_feature` removed). Metadata in a different order was reordered correctly.
- **Normalized data check:** a TPM table (non-integer values) and a log-transformed table (negative values) were refused with an explanation. Estimated counts (RSEM/Salmon-style decimals) were accepted and rounded only when **Round estimated counts** was ticked.
- **Sample sheet editor:** counts uploaded without metadata opened the Sample Sheet tab. A `Condition` column was added and filled in inside the app, and the analysis then ran normally.
- **Load from GEO:** tested against local copies of NCBI's file formats (raw counts `.tsv.gz` with Entrez IDs, Series Matrix files, the matrix folder listing, a bot-check page and missing files). The sample sheet was built from the sample characteristics (`cell_line`, `treatment`). A sample without counts was left out with a note. Organism and ID type were set automatically. A mouse study and a study split by platform (`-GPL…_series_matrix`) loaded correctly. Missing counts and bot-check pages gave a clear message. The analysis then ran end to end, with HYPOXIA again the top pathway.
