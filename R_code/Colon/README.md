# Colon Mix — Spatial Transcriptomics: Bacterial Microbiome Profiling

Spatial transcriptomics analysis of a defined bacterial community (22 species) colonizing mouse colon, profiled by DBiT-seq. The pipeline maps per-species bacterial UMI distributions across tissue, compares bacterial abundance between the Mix and Germ-Free conditions, and validates spatial abundances against bulk metagenomics.

---

## Data

| File | Description |
|------|-------------|
| `GSM8921749_colonmix_stdata.tsv` | Spot × gene expression matrix (GEO accession GSM8921749) |
| `bulk_abundance.tsv` | Metagenomics relative abundance and mapped read counts per species |
| `../Colon_GF/prefix_totals.tsv` | Germ-Free condition prefix totals (produced by the companion Colon_GF analysis) |

---

## Requirements

**R ≥ 4.1**

Install required packages:

```r
# Core (required)
install.packages(c("ggplot2", "dplyr", "tidyr", "grid", "scales", "viridis"))

# Recommended (used if available)
install.packages(c("ggrepel", "pheatmap", "ggforce", "reshape2", "OpenImageR"))
```

---

## Usage

1. Set the working directory to the folder containing `GSM8921749_colonmix_stdata.tsv`.
2. Open `Colonmix_analysis.R` and edit **Section 0 (Configuration)** if needed — particularly `gf_prefix_totals` and `metagenomics_file` paths.
3. Source the script:

```r
source("Colonmix_analysis.R")
```

All outputs are written to `results/`.

---

## Outputs

```
results/
├── prefix_totals.tsv                  Total UMI per bacterial species
│
├── species_maps/
│   └── <PREFIX>_UMI.pdf               Spatial UMI heatmap for each species (22 files)
│
├── violin_UMI.pdf                     UMI distribution per species across all spots
├── UMI_per_species_barplot.pdf/.png   Total UMI bar chart (all species)
├── total_UMI_pie.pdf                  Pie chart of species proportions
├── colocalization_heatmap.pdf         Pearson co-localization heatmap (requires pheatmap)
├── spatial_pie_chart.pdf              Per-spot composition pie chart (requires ggforce)
│
├── Mix_vs_GF_scatter.pdf              Mix vs Germ-Free scatter plot (Spearman ρ)
├── compare_mix_vs_GF.tsv              Mix vs GF comparison table
│
├── metagenomics_vs_mix_scatter.pdf    Metagenomics vs spatial relative abundance scatter
└── metagenomics_vs_mix.tsv            Merged metagenomics + spatial abundance table
```

---

## Bacterial Species Panel (22 species)

| Prefix | Species |
|--------|---------|
| GOZ73 | *Akkermansia muciniphila* |
| CL2 | *Anaerostipes hadrus* |
| F2Y25 | *Bacteroides caccae* |
| K6V30 | *Bacteroides faecis* |
| I6J55 | *Bacteroides fragilis* |
| Bovatus | *Bacteroides ovatus* |
| HMPREF2141 | *Bacteroides uniformis* |
| RY69 | *Bifidobacterium breve* |
| DW767 | *Blautia obeum* |
| LA357 | *Clostridium sporogenes* |
| CSV91 | *Collinsella aerofaciens* |
| DWX94 | *Coprococcus eutactus* |
| AQ990 | *Dorea longicatena* |
| LUA64 | *Eggerthella lenta* |
| ENT | *Enterococcus faecalis* |
| ELI | *Eubacterium callanderi* |
| Fi14EGH31 | *Faecalibacillus intestinalis* |
| C4Q21 | *Faecalibacterium prausnitzii* |
| DXB85 | *Parabacteroides merdae* |
| FYC25 | *Paraprevotella clara* |
| DWX93 | *Roseburia hominis* |
| CDL20 | *Ruminococcus gnavus* |
