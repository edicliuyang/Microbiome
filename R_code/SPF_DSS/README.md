# SPF vs IBD Spatial Transcriptomics Analysis

R scripts for analyzing host-microbiome interactions in spatial transcriptomics data comparing Specific Pathogen Free (SPF) and Inflammatory Bowel Disease (IBD) conditions.

## Data

Raw spatial transcriptomics data (ST-seq) is available from GEO:
- `GSM8921752` — SPF sample (`GSM8921752_SPF_stdata.tsv`)
- `GSM8921753` — IBD sample (`GSM8921753_SPFIBD_stdata.tsv`)

Download both files and place them in the working directory before running the scripts.

## Scripts

Run scripts in order:

| Script | Description |
|--------|-------------|
| `1_click_new.R` | Load and preprocess SPF and IBD spatial data; identify bacterial gene prefixes |
| `2_DEG_SPF_vs_IBD.R` | Differential expression analysis between SPF and IBD spots |
| `3_spatial_custom_filtered.R` | Spatial visualization of selected gene sets |
| `4_host_bacteria_colocalization.R` | Host-bacteria colocalization analysis |
| `5_host_bacteria_additional_figures.R` | Additional colocalization figures |
| `6_publication_figures.R` | Publication-quality figure generation |

## Dependencies

Install required R packages:

```r
install.packages(c("Seurat", "ggplot2", "dplyr", "ggrepel", "cowplot",
                   "gridExtra", "viridis", "ggpubr", "Matrix", "tidyr"))
```
