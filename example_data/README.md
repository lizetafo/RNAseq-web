# Example data

Two ready-to-use datasets. Each folder has a raw count matrix (`counts.csv`, genes × samples) and a sample sheet (`metadata.csv`).

| Folder | Organism | Gene IDs | Samples | Comparison | App settings |
|---|---|---|---|---|---|
| `prostate_hypoxia/` | Human | Ensembl (`ENSG…`) | 8 (4 + 4) | PC3 and LNCaP prostate cancer cells, **Hypoxia vs Normoxia** | Human · ENSEMBL · design `~Condition` · control group `Normoxia` · contrast `Condition, Hypoxia, Normoxia`. Also try `~CellLine + Condition` (groups `Condition, CellLine`) |
| `tumor_vs_normal/` | Human | Ensembl | 37 (18 cancer + 19 normal) | **colorectal cancer vs normal colon** | Human · ENSEMBL · `~type` · control group `normal` · contrast `type, cancer, normal` |

`prostate_hypoxia/metadata.csv` includes a `CellLine` column (LNCaP or PC3, taken from the sample names), so you can compare a one-factor design (`~Condition`) with a two-factor design that accounts for cell line (`~CellLine + Condition`). The **Quality Control** tab shows why the second is better for these data.

## Data sources
Public datasets from NCBI GEO, reduced to the samples used here:

| Folder | GEO series | Samples used |
|---|---|---|
| `prostate_hypoxia/` | [GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305): Guo et al., *Nature Communications* (2019) | LNCaP (empty vector) and PC3 (siCtrl), normoxia and hypoxia, 2 replicates each |
| `tumor_vs_normal/` | [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760): colorectal cancer | 18 primary colorectal tumors (SRR975551–SRR975568) and 19 normal colon samples (SRR975569–SRR975587) |

## Input format

**Counts:** first column = gene IDs; one column per sample; raw integer counts (for example from featureCounts or HTSeq). Do not use normalized values such as TPM or FPKM.
```
GeneID,Sample1,Sample2,Sample3
ENSG00000000003,348,609,693
```

**Metadata:** first column = sample names, matching the count column names **exactly and in the same order**; then one column per experimental factor.
```
ID,Condition
Sample1,Normoxia
Sample2,Hypoxia
```
CSV, TSV and TXT files are all accepted.
