# Example data

Two ready-to-use datasets. Each folder has a raw count matrix (`counts.csv`, genes × samples) and a sample sheet (`metadata.csv`).

| Folder | Organism | Gene IDs | Samples | Comparison | App settings |
|---|---|---|---|---|---|
| `prostate_hypoxia/` | Human | Ensembl (`ENSG…`) | 8 (4 + 4) | PC3 and LNCaP prostate cancer cells, **Hypoxia vs Normoxia** | Human · ENSEMBL · design `~Condition` · control group `Normoxia` · contrast `Condition, Hypoxia, Normoxia`. Also try `~CellLine + Condition` (groups `Condition, CellLine`) |
| `tumor_vs_normal/` | Human | Ensembl | 37 (18 cancer + 19 normal) | **colorectal cancer vs normal colon** | Human · ENSEMBL · `~type` · control group `normal` · contrast `type, cancer, normal` |

`prostate_hypoxia/metadata.csv` includes a `CellLine` column (LNCaP or PC3, taken from the sample names), so you can compare a one-factor design (`~Condition`) with a two-factor design that accounts for cell line (`~CellLine + Condition`). The **Quality Control** tab shows why the second is better for these data. To analyze one cell line on its own, as the original study did for PC3, use **Analyze a subset of samples** → `CellLine` → `PC3`. Note that each group has only 2 replicates, so the app shows a low-replication warning.

## Data sources
Public datasets from NCBI GEO, reduced to the samples used here:

| Folder | GEO series | Samples used |
|---|---|---|
| `prostate_hypoxia/` | [GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305): Guo et al., *Nature Communications* (2019) | LNCaP (empty vector) and PC3 (siCtrl), normoxia and hypoxia, 2 replicates each |
| `tumor_vs_normal/` | [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760): colorectal cancer | 18 primary colorectal tumors (SRR975551–SRR975568) and 19 normal colon samples (SRR975569–SRR975587) |

**Please cite the original studies if you use these data:**
- Guo H, Ci X, Ahmed M, et al. ONECUT2 is a driver of neuroendocrine prostate cancer. *Nature Communications*. 2019;10(1):278. https://doi.org/10.1038/s41467-018-08133-6
- Kim SK, Kim SY, Kim JH, et al. A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. *Molecular Oncology*. 2014;8(8):1653–1666. https://doi.org/10.1016/j.molonc.2014.06.016

## Input format

**Counts:** first column = gene IDs; one column per sample; raw integer counts (for example from featureCounts or HTSeq). Normalized values such as TPM or FPKM are refused, because DESeq2 needs raw counts.

The app also accepts featureCounts output, several per-sample HTSeq or STAR count files (merged automatically), a GEO `_RAW.tar` bundle, and `.gz` files. You can also skip files entirely and use **Load from GEO**.
```
GeneID,Sample1,Sample2,Sample3
ENSG00000000003,348,609,693
```

**Metadata (optional):** first column = sample names, matching the count column names **exactly** (a different order is fixed automatically); then one column per experimental factor. Without a metadata file, fill in the sample sheet in the app's **Sample Sheet** tab.
```
ID,Condition
Sample1,Normoxia
Sample2,Hypoxia
```
CSV, TSV and TXT files are all accepted.
