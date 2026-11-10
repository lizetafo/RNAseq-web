# User guide

## 1. One-time setup
Open RStudio, set the working directory to this repository folder, and run:
```r
source("setup/install_packages.R")
```
Then go to **Session → Restart R**.

> Organism databases are large downloads. The setup script installs the **human** database. Others (mouse, rat, fly, zebrafish, yeast) are installed automatically the first time you select them in the app, which can take a few minutes once.

## 2. Start the app
Open `app.R` in RStudio and click **Run App**.

---

## 3. Upload your data and settings

![Upload and settings](screenshots/1_upload_and_settings.png)

1. **Upload Counts File**: raw counts, genes × samples (CSV, TSV or TXT). See [input format](../example_data/README.md#input-format).
2. **Upload Metadata File**: one row per sample, in the same order as the count columns.
3. Leave **Convert NA in counts to zero** ticked unless you want the app to stop on missing values.
   **No data?** Click **Or load example data (prostate hypoxia)** instead of steps 1–4.

4. Click **Upload**. The app checks that sample names match, removes stray spaces from group labels, then detects the grouping column (the first metadata column with 2–10 groups) and fills in the settings below.

5. **Choose the control group** in **Reference level (control group)**. The app does not guess it. The **contrast** then fills in automatically (case first, control last). Check the other settings:

| Setting | What it means | Example |
|---|---|---|
| **Organism Database** | Used to add gene names | Human (org.Hs.eg.db) |
| **ID Type** | The type of gene IDs in your counts file: `ENSG…` = ENSEMBL, numbers = ENTREZID, gene names = SYMBOL | ENSEMBL |
| **Experimental Grouping Factor(s)** | Metadata column(s) that define groups; separate several with commas | `Condition` |
| **Reference Column / Reference level (control group)** | The **control** group that fold changes are measured against. Choose it from the dropdown. | `Condition` / `Normoxia` |
| **Experimental Design Formula** | DESeq2 design; add covariates with `+` | `~Condition` or `~Activation + Infection` |
| **Contrast (Factor, Case, Control)** | The comparison to report: **case first, control last** | `Condition, Hypoxia, Normoxia` |

With contrast `Condition, Hypoxia, Normoxia`, a **positive log2 fold change means higher in Hypoxia**.

**Thresholds**
- **Adjusted p-value threshold (alpha)**: significance cut-off used by DESeq2 (default 0.05).
- **log2FC threshold** and **adjusted p-value threshold** for the plots: which genes are colored as significant on the volcano plot.
- **Transformation for heatmap**: normalized counts, VST or rlog. VST and rlog stabilize the variance of low counts, and VST is a good default.
- **Number of top genes**: how many genes to label and show in the heatmap and bar plot.

6. Click **Run**. If no control group is chosen, the app asks for one. DESeq2 usually takes from a few seconds to about a minute.

> **Privacy:** the online version runs on Posit Connect Cloud, so uploaded files are processed on Posit's servers during your session. The app does not save them. Do not upload confidential, unpublished or patient-identifiable data online; run the app on your own computer instead.

---

## 4. Quality control (check this first)

The **Quality Control** tab shows whether the samples and the model look healthy **before** you interpret results. Use **Color / group samples by** to color samples by any metadata column.

![PCA](screenshots/6_qc_pca.png)

| Plot | What to look for |
|---|---|
| **PCA** | Replicates of the same group should cluster together. If samples separate mainly by something else (cell line, batch, donor), add that column to the design, for example `~CellLine + Condition`. |
| **Sample Distances** | Darker squares = more similar samples. Samples should group by condition. A sample unlike all others may be an outlier. |
| **Library Size** | Total reads per sample and the DESeq2 size factor (SF). Big differences are normal, because DESeq2 corrects for them, but a sample with very few reads deserves a second look. |
| **Dispersion** | DESeq2's variability estimates. Blue (final) points should follow the red trend line. |
| **P-value Histogram** | Ideally flat with a peak near 0. A histogram rising toward 1 suggests a source of variation missing from the design. |

**Example:** in the prostate dataset, the PCA shows that **cell line** (LNCaP vs PC3) explains 98.6% of the variation. With `~Condition`, 302 genes are significant. With the design `~CellLine + Condition` (Experimental Grouping Factor(s): `Condition, CellLine`), 4,895 genes are significant, because the cell line difference no longer hides the hypoxia effect.

---

## 5. Result plots

**Volcano plot:** every gene, log2 fold change vs. significance, with the top genes labeled. A table of the top 15 genes appears beside it.

![Volcano plot](screenshots/2_volcano_plot.png)

**Heatmap:** the top genes across all samples, scaled per gene (z-score) and clustered.

![Heatmap](screenshots/3_heatmap.png)

**Bar plot:** the most significant genes (lowest adjusted p-value) shown as bars of their log2 fold change, colored up (pink) or down (blue).

![Bar plot](screenshots/4_bar_plot.png)

**MA plot:** fold change against average expression for every gene. Genes with low counts tend to show larger, less reliable fold changes.

![MA plot](screenshots/8_ma_plot.png)

Each plot has a **Download** button.

## 6. Gene Explorer

Pick any gene, by typing its name or ID, to see its normalized counts in each group and its DESeq2 statistics. The most significant gene is shown first.

![Gene explorer](screenshots/9_gene_explorer.png)

## 7. Results tables

![Results tables](screenshots/5_results_tables.png)

| Tab | Contents | Download button |
|---|---|---|
| **Annotated** | DESeq2 results + gene names | Download Annotated Results |
| **Raw DESeq2 Results** | DESeq2 results with gene IDs only | Download DESeq2 Results |
| **Normalized Counts** | DESeq2 size-factor-normalized counts | Download Normalized Counts |

Column meanings are in [outputs.md](outputs.md).

---

## Troubleshooting

| Message or problem | Fix |
|---|---|
| *Some samples in counts data missing from metadata* | A count column name has no matching row in the metadata. Check spelling and spaces. |
| *Mismatch in sample order between counts and metadata* | Sort the metadata rows in the same order as the count columns. |
| *NA values in counts* | Tick "Convert NA in counts to zero", or fix the file. |
| *Design formula is not valid* | Start the formula with `~`, for example `~Condition`. |
| *Invalid contrast entered* | Use exactly three items: `factor, case, control`. |
| More groups than expected | Look for typos in group labels, for example `Activated` vs `activated`. Stray spaces are removed automatically. |
| Few genes get gene names | The ID type or organism is probably wrong. |
| Package installation fails on first run | Run `setup/install_packages.R` first, then restart R. |
| *Please choose the reference (control) group* | Pick the control group in the **Reference level** dropdown, then click Run. |
| A message says an organism database is installing | One-time download the first time you use that organism. Wait for it to finish. |
