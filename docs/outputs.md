# Outputs

## Results tables (Annotated / Raw DESeq2 Results)
One row per gene that passed the low-count filter (total count across all samples ≥ 10).

| Column | Meaning |
|---|---|
| *(row name)* | Gene ID from your counts file |
| `baseMean` | Mean of the normalized counts across all samples |
| `log2FoldChange` | log2 (case ÷ control), for the contrast you entered. **Positive = higher in the case group.** |
| `lfcSE` | Standard error of the log2 fold change |
| `stat` | Wald test statistic |
| `pvalue` | Raw p-value |
| `padj` | Benjamini–Hochberg adjusted p-value (false discovery rate). **Use this for significance.** |
| `GeneName` | Gene symbol from the organism database *(Annotated table only)*. Empty when the ID has no symbol, which is common for non-coding genes. |

`padj` is `NA` when DESeq2 excludes a gene from testing: all counts are zero, it has an extreme outlier, or its mean count is too low for independent filtering. This is normal DESeq2 behavior.

## Normalized counts
Counts divided by DESeq2 size factors, so samples sequenced to different depths can be compared. Use them for plots of individual genes. Don't feed them back into DESeq2, which always needs raw counts.

## Plots
| Plot | What it shows |
|---|---|
| **Volcano** | x = log2 fold change, y = −log10 adjusted p-value. Dashed lines show the log2FC and adjusted p-value thresholds. The top genes are labeled. |
| **Heatmap** | Top genes by adjusted p-value × samples, scaled per gene (z-score), with genes and samples clustered. Choose normalized counts, VST or rlog values. |
| **Bar plot** | Log2 fold change of the top genes by adjusted p-value, colored by direction |
| **MA plot** | x = log10 mean normalized count, y = log2 fold change; colored like the volcano plot |

## Quality-control plots
| Plot | What it shows |
|---|---|
| **PCA** | Principal component analysis of the VST-transformed data, using the most variable genes (default 500, adjustable). Axis labels give the % of variance explained. |
| **Sample distances** | Euclidean distances between samples on VST data, clustered, with samples annotated by the chosen metadata column |
| **Library size** | Total raw counts per sample (millions), labeled with the DESeq2 size factor |
| **Dispersion** | `DESeq2::plotDispEsts`: gene-wise (black), fitted (red) and final (blue) dispersion estimates |
| **P-value histogram** | Distribution of raw p-values (bins of 0.025) |

## Gene Explorer
Normalized counts (+0.5, log scale) of one chosen gene per sample, grouped by the first experimental factor, with its DESeq2 statistics (baseMean, log2FoldChange, lfcSE, pvalue, padj).
