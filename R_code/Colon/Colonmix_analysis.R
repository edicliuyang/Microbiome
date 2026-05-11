################################################################################
# Colon Mix — Spatial Transcriptomics: Bacterial Microbiome Profiling
#
# Required input files (place in working directory):
#   GSM8921749_colonmix_stdata.tsv   spot × gene expression matrix
#   ../Colon_GF/prefix_totals.tsv    GF condition prefix totals (optional)
#   bulk_abundance.tsv               metagenomics read counts    (optional)
################################################################################


## ============================================================
## 0) CONFIGURATION — edit before running
## ============================================================

cfg <- list(

  stdata_file      = "GSM8921749_colonmix_stdata.tsv",
  gf_prefix_totals = "../Colon_GF/prefix_totals.tsv",  # NULL to skip
  metagenomics_file = "bulk_abundance.tsv",             # NULL to skip
  image_file       = NULL,   # optional tissue background image

  xlim = c(0, 51),
  ylim = c(51, 1),
  point_size = 3.5,

  # Bacterial species prefixes (must match gene columns: <PREFIX>_RS...)
  prefixes = c(
    "GOZ73", "CL2",        "F2Y25",    "K6V30",      "I6J55",     "Bovatus",
    "HMPREF2141", "RY69",  "DW767",    "LA357",      "CSV91",     "DWX94",
    "AQ990",  "LUA64",     "ENT",      "ELI",        "Fi14EGH31", "C4Q21",
    "DXB85",  "FYC25",     "DWX93",    "CDL20"
  ),

  # Abbreviated labels for plots
  species_labels = c(
    "GOZ73" = "A. muciniphila", "CL2" = "A. hadrus",
    "F2Y25" = "B. caccae",     "K6V30" = "B. faecis",
    "I6J55" = "B. fragilis",   "Bovatus" = "B. ovatus",
    "HMPREF2141" = "B. uniformis", "RY69" = "B. breve",
    "DW767" = "B. obeum",      "LA357" = "C. sporogenes",
    "CSV91" = "C. aerofaciens","DWX94" = "C. eutactus",
    "AQ990" = "D. longicatena","LUA64" = "E. lenta",
    "ENT"   = "E. faecalis",   "ELI"   = "E. callanderi",
    "Fi14EGH31" = "F. intestinalis", "C4Q21" = "F. prausnitzii",
    "DXB85" = "P. merdae",     "FYC25" = "P. clara",
    "DWX93" = "R. hominis",    "CDL20" = "R. gnavus"
  ),

  # Full names for merging with metagenomics table
  species_full = c(
    "GOZ73" = "Akkermansia_muciniphila",  "CL2" = "Anaerostipes_hadrus",
    "F2Y25" = "Bacteroides_caccae",       "K6V30" = "Bacteroides_faecis",
    "I6J55" = "Bacteroides_fragilis",     "Bovatus" = "Bacteroides_ovatus",
    "HMPREF2141" = "Bacteroides_uniformis","RY69" = "Bifidobacterium_breve",
    "DW767" = "Blautia_obeum",            "LA357" = "Clostridium_sporogenes",
    "CSV91" = "Collinsella_aerofaciens",  "DWX94" = "Coprococcus_eutactus",
    "AQ990" = "Dorea_longicatena",        "LUA64" = "Eggerthella_lenta",
    "ENT"   = "Enterococcus_faecalis",    "ELI"   = "Eubacterium_callanderi",
    "Fi14EGH31" = "Faecalibacillus_intestinalis",
    "C4Q21" = "Faecalibacterium_prausnitzii",
    "DXB85" = "Parabacteroides_merdae",   "FYC25" = "Paraprevotella_clara",
    "DWX93" = "Roseburia_hominis",        "CDL20" = "Ruminococcus_gnavus"
  ),

  violin_ylim  = 300,
  cor_threshold = 0.3,

  # Seurat PCA settings (for elbow plot)
  seurat_dims = 1:10,   # adjust after inspecting the elbow plot
  out_dir = "results"
)


## ============================================================
## 1) PACKAGES
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(scales)
  library(viridis)
})

has_pkg <- function(x) requireNamespace(x, quietly = TRUE)
if (has_pkg("OpenImageR")) library(OpenImageR)
if (has_pkg("ggforce"))    library(ggforce)
if (has_pkg("pheatmap"))   library(pheatmap)
if (has_pkg("ggrepel"))    library(ggrepel)
if (has_pkg("reshape2"))   library(reshape2)


## ============================================================
## 2) HELPERS
## ============================================================

msg       <- function(...) cat(sprintf(...), "\n")
safe_mkdir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
`%||%`    <- function(a, b) if (is.null(a)) b else a

get_bg_grob <- function(path) {
  if (is.null(path) || !file.exists(path) || !has_pkg("OpenImageR")) return(NULL)
  rasterGrob(OpenImageR::readImage(path),
             width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
}

add_bg <- function(p, grob) {
  if (is.null(grob)) return(p)
  p + annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

void_spatial <- function() {
  theme_void() +
  theme(plot.title   = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))
}

robust_upper <- function(x) {
  x <- x[x > 0 & is.finite(x)]
  if (length(x) == 0) return(1)
  max(if (length(x) >= 10) quantile(x, 0.99) * 1.5 else max(x) * 1.5, 1)
}

# Sum UMI per spot per prefix (genes are named <PREFIX>_RS...)
sum_prefix_umi <- function(df, prefixes) {
  out <- df %>% dplyr::select(A, B)
  for (pref in prefixes) {
    cols <- grep(paste0("^", pref, "_"), colnames(df), value = TRUE)
    out[[pref]] <- if (length(cols) > 0) rowSums(df[, cols, drop = FALSE], na.rm = TRUE) else 0
  }
  out
}

# Long format with per-spot fractions + arc angles for pie charts
build_umi_long <- function(spu, prefixes) {
  spu %>%
    pivot_longer(all_of(prefixes), names_to = "Prefix", values_to = "UMI") %>%
    group_by(A, B) %>%
    mutate(Total_UMI = sum(UMI),
           Fraction  = if_else(Total_UMI > 0, UMI / Total_UMI, 0)) %>%
    arrange(A, B, Prefix) %>%
    mutate(start = 2 * pi * cumsum(lag(Fraction, default = 0)),
           end   = 2 * pi * cumsum(Fraction)) %>%
    ungroup()
}

default_colors <- function(n) {
  pal <- c("#F0CE58","#B487B7","#289E92","#EB545C","#5084C2","#DBA091",
           "#D7EF9B","#EF7512","#6A5ACD","#20B2AA","#FF69B4","#87CEEB",
           "#FFD700","#A52A2A","#40E0D0","#D2691E","#9ACD32","#00CED1",
           "#DC143C","#8A2BE2","#556B2F","#FF8C00")
  rep(pal, length.out = n)
}


## ============================================================
## 3) LOAD DATA
## ============================================================

stopifnot(file.exists(cfg$stdata_file))
safe_mkdir(cfg$out_dir)
bg <- get_bg_grob(cfg$image_file)

msg("Loading %s ...", cfg$stdata_file)
raw <- read.table(cfg$stdata_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
if (!("X" %in% colnames(raw))) colnames(raw)[1] <- "X"

# Parse spot coordinates
df <- raw %>%
  tidyr::separate(X, into = c("A", "B"), sep = "x",
                  remove = FALSE, convert = TRUE)

# Build per-spot prefix UMI matrix
spot_pumi  <- sum_prefix_umi(df, cfg$prefixes)
prefix_colors <- setNames(default_colors(length(cfg$prefixes)), cfg$prefixes)

# Global totals (used in sections 5 & 6)
totals_df <- spot_pumi %>%
  summarise(across(all_of(cfg$prefixes), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "Prefix", values_to = "Total_UMI") %>%
  mutate(Fraction = Total_UMI / sum(Total_UMI),
         Species  = cfg$species_labels[Prefix]) %>%
  arrange(desc(Fraction)) %>%
  mutate(Species = factor(Species, levels = Species))

write.table(totals_df, file.path(cfg$out_dir, "prefix_totals.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
msg("Loaded %d spots, %d bacterial genes detected.",
    nrow(df),
    sum(sapply(cfg$prefixes, function(p)
      sum(grepl(paste0("^", p, "_"), colnames(df))))))


## ============================================================
## 4) PCA ELBOW PLOT (host genes)
## ============================================================

# Host genes: exclude bacterial, mitochondrial, and Gm pseudogenes
gene_cols  <- setdiff(colnames(raw), "X")
bact_genes <- grep(paste0("^(", paste(cfg$prefixes, collapse = "|"), ")_"),
                   gene_cols, value = TRUE)
host_cols  <- setdiff(gene_cols, c(bact_genes, "unmapped"))
host_cols  <- host_cols[!grepl("^(Gm|mt-|Mt-|MT-|Rn18s|Rn28s|Rn5s)", host_cols)]

rownames(raw) <- raw$X

# Keep only spots with at least 1 host gene detected (avoids log(0) in SCTransform)
host_counts <- rowSums(raw[, host_cols, drop = FALSE])
keep_spots  <- names(host_counts)[host_counts > 0]
msg("%d / %d spots retained after removing zero-count spots.",
    length(keep_spots), nrow(raw))

mat  <- Matrix::Matrix(as.matrix(t(raw[keep_spots, host_cols])), sparse = TRUE)
pbmc <- CreateSeuratObject(mat, min.cells = 1, min.features = 1,
                           project = "ColonMix")
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)

# Variance explained by each PC
pca_stdev    <- pbmc[["pca"]]@stdev
pct_var      <- pca_stdev^2 / sum(pca_stdev^2) * 100
cum_var_10   <- round(sum(pct_var[1:max(cfg$seurat_dims)]), 1)
msg("First %d PCs explain %.1f%% of total variance.", max(cfg$seurat_dims), cum_var_10)

p_elbow <- ElbowPlot(pbmc, ndims = 30) +
  geom_vline(xintercept = max(cfg$seurat_dims), linetype = "dashed",
             color = "red", linewidth = 0.8) +
  annotate("text", x = max(cfg$seurat_dims) + 0.5, y = Inf,
           label = paste0("PC1–", max(cfg$seurat_dims), ": ",
                          cum_var_10, "% variance"),
           hjust = 0, vjust = 1.5, color = "red", size = 4) +
  labs(title = "PCA Elbow Plot — host genes",
       x = "Principal component", y = "Standard deviation") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(cfg$out_dir, "PCA_elbow_plot.pdf"),
       plot = p_elbow, width = 7, height = 5)
ggsave(file.path(cfg$out_dir, "PCA_elbow_plot.png"),
       plot = p_elbow, width = 7, height = 5, dpi = 300)

msg("Elbow plot saved. Inspect it and update cfg$seurat_dims if needed.")

# Clustering (uses cfg$seurat_dims set above)
pbmc <- RunUMAP(pbmc, dims = cfg$seurat_dims, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = cfg$seurat_dims, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = cfg$seurat_resolution, verbose = FALSE)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                               min.pct = 0.1, logfc.threshold = 0.25,
                               verbose = FALSE)
write.table(pbmc.markers, file.path(cfg$out_dir, "cluster_markers.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("Clustering done: %d clusters found.", length(levels(Idents(pbmc))))


## ============================================================
## 4b) GO & KEGG PATHWAY ENRICHMENT (host cluster markers)
##
## Requires:
##   BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))
##   Use org.Hs.eg.db instead if analysing human data.
## ============================================================

if (has_pkg("clusterProfiler") && has_pkg("org.Mm.eg.db")) {

  library(clusterProfiler)
  library(org.Mm.eg.db)

  safe_mkdir(file.path(cfg$out_dir, "enrichment"))

  # Significant markers across all clusters
  sig_markers <- pbmc.markers %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    filter(!grepl("^(Gm|mt-|Rn18s)", gene))

  gene_df <- clusterProfiler::bitr(
    unique(sig_markers$gene),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )

  # ---- GO enrichment (Biological Process) ----
  go_res <- clusterProfiler::enrichGO(
    gene          = gene_df$ENTREZID,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  # ---- KEGG enrichment ----
  kegg_res <- clusterProfiler::enrichKEGG(
    gene          = gene_df$ENTREZID,
    organism      = "mmu",        # mmu = mouse; hsa = human
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  kegg_res <- clusterProfiler::setReadable(kegg_res, OrgDb = org.Mm.eg.db,
                                           keyType = "ENTREZID")

  # ---- Save result tables ----
  write.table(as.data.frame(go_res),
              file.path(cfg$out_dir, "enrichment", "GO_BP_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(as.data.frame(kegg_res),
              file.path(cfg$out_dir, "enrichment", "KEGG_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Visualize: dotplots ----
  n_show <- 20   # number of terms to display

  p_go_dot <- dotplot(go_res, showCategory = n_show) +
    ggtitle("GO Enrichment — Biological Process") +
    theme_bw(base_size = 13) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))

  p_kegg_dot <- dotplot(kegg_res, showCategory = n_show) +
    ggtitle("KEGG Pathway Enrichment") +
    theme_bw(base_size = 13) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))

  # ---- Visualize: barplots (ggplot2) ----
  go_top <- as.data.frame(go_res) %>%
    slice_head(n = n_show) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  p_go_bar <- ggplot(go_top, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#d73027", high = "#91bfdb",
                        name = "Adjusted\np-value") +
    labs(title = "GO Enrichment — Biological Process",
         x = "Gene count", y = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))

  kegg_top <- as.data.frame(kegg_res) %>%
    slice_head(n = n_show) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  p_kegg_bar <- ggplot(kegg_top, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#d73027", high = "#91bfdb",
                        name = "Adjusted\np-value") +
    labs(title = "KEGG Pathway Enrichment",
         x = "Gene count", y = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))

  # ---- Save plots ----
  ggsave(file.path(cfg$out_dir, "enrichment", "GO_dotplot.pdf"),
         plot = p_go_dot,   width = 10, height = 8)
  ggsave(file.path(cfg$out_dir, "enrichment", "GO_barplot.pdf"),
         plot = p_go_bar,   width = 10, height = 8)
  ggsave(file.path(cfg$out_dir, "enrichment", "KEGG_dotplot.pdf"),
         plot = p_kegg_dot, width = 10, height = 8)
  ggsave(file.path(cfg$out_dir, "enrichment", "KEGG_barplot.pdf"),
         plot = p_kegg_bar, width = 10, height = 8)

  msg("GO and KEGG enrichment saved to results/enrichment/")

} else {
  msg("clusterProfiler or org.Mm.eg.db not installed — skipping enrichment.")
  msg("Install with: BiocManager::install(c('clusterProfiler','org.Mm.eg.db'))")
}


## ============================================================
## 5) PER-SPECIES SPATIAL UMI MAPS  (update cfg$seurat_dims first if needed)
## ============================================================

safe_mkdir(file.path(cfg$out_dir, "species_maps"))

for (pref in cfg$prefixes) {
  umi_vec <- spot_pumi[[pref]]
  if (all(umi_vec == 0)) { msg("No UMI for %s — skip.", pref); next }

  dfp   <- spot_pumi %>% dplyr::select(A, B) %>% mutate(UMI = umi_vec)
  upper <- robust_upper(umi_vec)
  label <- cfg$species_labels[[pref]]

  p <- ggplot(dfp, aes(x = as.numeric(A), y = as.numeric(B), color = UMI)) +
    geom_point(shape = 16, size = cfg$point_size) +
    scale_color_gradient(low = "white", high = "red",
                         limits = c(0, upper), oob = scales::squish) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 28)) +
    ggtitle(label) +
    coord_equal(xlim = cfg$xlim, ylim = cfg$ylim) +
    scale_y_reverse() +
    void_spatial()

  ggsave(file.path(cfg$out_dir, "species_maps", paste0(pref, "_UMI.pdf")),
         plot = add_bg(p, bg), width = 8.6, height = 8.6)
}


## ============================================================
## 5) COMPOSITION PLOTS
## ============================================================

umi_long <- build_umi_long(spot_pumi, cfg$prefixes)

# --- Violin: UMI distribution per species ---
p_vln <- ggplot(umi_long, aes(x = Prefix, y = UMI, fill = Prefix)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  scale_fill_manual(values = prefix_colors) +
  scale_x_discrete(labels = cfg$species_labels[cfg$prefixes]) +
  scale_y_continuous(limits = c(0, cfg$violin_ylim),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "UMI per spot",
       title = "Bacterial species UMI distribution") +
  theme_classic(base_size = 16) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1,
                                     size = 11, face = "bold.italic"),
        plot.title    = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

ggsave(file.path(cfg$out_dir, "violin_UMI.pdf"),
       plot = p_vln, width = 12, height = 6)

# --- Bar chart: total UMI per species ---
p_bar <- ggplot(totals_df, aes(x = Species, y = Total_UMI)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = Total_UMI), vjust = -0.3, size = 3.5) +
  labs(title = "Total UMIs per species", x = NULL, y = "Total UMI") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

ggsave(file.path(cfg$out_dir, "UMI_per_species_barplot.pdf"),
       plot = p_bar, width = 10, height = 6)
ggsave(file.path(cfg$out_dir, "UMI_per_species_barplot.png"),
       plot = p_bar, width = 10, height = 6, dpi = 300)

# --- Global pie chart ---
p_pie <- ggplot(totals_df, aes(x = "", y = Fraction, fill = Prefix)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = prefix_colors,
                    labels = cfg$species_labels[cfg$prefixes]) +
  geom_text(aes(label = ifelse(Fraction > 0.02,
                               paste0(round(Fraction * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 4) +
  theme_void() +
  ggtitle("Total UMI proportion per species") +
  theme(plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_blank())

ggsave(file.path(cfg$out_dir, "total_UMI_pie.pdf"),
       plot = p_pie, width = 7, height = 7)

# --- Co-localization heatmap (requires pheatmap) ---
if (has_pkg("pheatmap")) {
  mat_c <- spot_pumi %>% dplyr::select(all_of(cfg$prefixes)) %>% as.matrix()
  colnames(mat_c) <- cfg$species_labels[colnames(mat_c)]
  cor_m <- cor(mat_c, method = "pearson", use = "pairwise.complete.obs")
  cor_m[is.na(cor_m)] <- 0

  pdf(file.path(cfg$out_dir, "colocalization_heatmap.pdf"), width = 9, height = 8)
  pheatmap::pheatmap(
    cor_m,
    cluster_rows = TRUE, cluster_cols = TRUE,
    display_numbers = TRUE, number_format = "%.2f",
    color  = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101),
    fontsize = 11, border_color = NA,
    main = "Co-localization heatmap (Pearson r)"
  )
  dev.off()
}

# --- Spatial pie chart (requires ggforce) ---
if (has_pkg("ggforce")) {
  umi_long_f <- umi_long %>%
    filter(Total_UMI > 5) %>%
    mutate(radius = sqrt(Total_UMI) * 0.04)

  p_spie <- ggplot(umi_long_f) +
    ggforce::geom_arc_bar(aes(
      x0 = as.numeric(A), y0 = as.numeric(B),
      r0 = 0, r = radius, start = start, end = end, fill = Prefix
    ), color = NA, alpha = 0.6) +
    scale_fill_manual(values = prefix_colors,
                      labels = cfg$species_labels[cfg$prefixes]) +
    coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
    scale_y_reverse() +
    void_spatial() +
    ggtitle("Bacterial composition per spot")

  ggsave(file.path(cfg$out_dir, "spatial_pie_chart.pdf"),
         plot = add_bg(p_spie, bg), width = 8.6, height = 8.6)
}


## ============================================================
## 6) MIX vs GF — SCATTER PLOT
## ============================================================

if (!is.null(cfg$gf_prefix_totals) && file.exists(cfg$gf_prefix_totals)) {

  gf_raw <- read.table(cfg$gf_prefix_totals, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
  colnames(gf_raw)[1:2] <- c("prefix", "GF_total_UMI")

  compare_df <- totals_df %>%
    dplyr::select(prefix = Prefix, mix_total_UMI = Total_UMI) %>%
    merge(gf_raw[, c("prefix","GF_total_UMI")], by = "prefix", all = TRUE)
  compare_df[is.na(compare_df)] <- 0
  compare_df$species_label <- cfg$species_labels[compare_df$prefix]

  # Spearman correlation
  cr <- cor.test(log10(compare_df$GF_total_UMI  + 1),
                 log10(compare_df$mix_total_UMI + 1),
                 method = "spearman")
  cor_text <- paste0("Spearman ρ = ", round(cr$estimate, 3),
                     "\nP = ", signif(cr$p.value, 3))

  p_gf <- ggplot(compare_df,
                 aes(x = log10(GF_total_UMI + 1),
                     y = log10(mix_total_UMI + 1))) +
    geom_point(size = 3, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
    {if (has_pkg("ggrepel"))
       ggrepel::geom_text_repel(aes(label = species_label),
                                size = 4, fontface = "italic")
     else
       geom_text(aes(label = species_label), size = 3,
                 vjust = -0.6, fontface = "italic")} +
    annotate("text",
             x = min(log10(compare_df$GF_total_UMI  + 1)) + 0.1,
             y = max(log10(compare_df$mix_total_UMI + 1)) - 0.1,
             label = cor_text, hjust = 0, size = 5) +
    theme_bw(base_size = 14) +
    labs(title = "Bacterial abundance: Mix vs Germ-Free",
         x = "Germ-Free total UMI (log₁₀ + 1)",
         y = "Mix total UMI (log₁₀ + 1)")

  ggsave(file.path(cfg$out_dir, "Mix_vs_GF_scatter.pdf"),
         plot = p_gf, width = 8, height = 7)

  write.table(compare_df,
              file.path(cfg$out_dir, "compare_mix_vs_GF.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("Mix vs GF scatter saved.")

} else {
  msg("GF prefix totals not found (%s) — skipping.",
      cfg$gf_prefix_totals %||% "NULL")
}


## ============================================================
## 7) METAGENOMICS vs MIX SPATIAL — SCATTER PLOT
## ============================================================

if (!is.null(cfg$metagenomics_file) && file.exists(cfg$metagenomics_file)) {

  meta <- read.table(cfg$metagenomics_file, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)

  mix_rel <- totals_df %>%
    mutate(species = cfg$species_full[Prefix],
           mix_rel_abund = Total_UMI / sum(Total_UMI)) %>%
    filter(!is.na(species)) %>%
    dplyr::select(species, mix_total_UMI = Total_UMI, mix_rel_abund)

  comp <- merge(meta, mix_rel, by = "species", all = TRUE)
  comp[is.na(comp)] <- 0
  comp$species_label <- sub("^([A-Za-z])[a-z]+_([A-Za-z]+)$", "\\1. \\2",
                             comp$species)
  comp <- comp[comp$relative_abundance > 0 | comp$mix_rel_abund > 0, ]

  comp$log10_meta <- log10(comp$relative_abundance + 1e-8)
  comp$log10_mix  <- log10(comp$mix_rel_abund       + 1e-8)

  cr2 <- cor.test(comp$log10_meta, comp$log10_mix, method = "spearman")
  cor_text2 <- paste0("Spearman ρ = ", round(cr2$estimate, 3),
                      "\nP = ", signif(cr2$p.value, 3))

  p_meta <- ggplot(comp, aes(x = log10_meta, y = log10_mix)) +
    geom_point(size = 3, color = "firebrick2") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
    {if (has_pkg("ggrepel"))
       ggrepel::geom_text_repel(aes(label = species_label),
                                size = 4, fontface = "italic")
     else
       geom_text(aes(label = species_label), size = 3,
                 vjust = -0.6, fontface = "italic")} +
    annotate("text",
             x = min(comp$log10_meta) + 0.15,
             y = max(comp$log10_mix)  - 0.15,
             label = cor_text2, hjust = 0, size = 5) +
    theme_bw(base_size = 14) +
    labs(title = "Relative abundance: metagenomics vs Mix spatial",
         x = "Metagenomics relative abundance (log₁₀)",
         y = "Mix spatial relative abundance (log₁₀)")

  ggsave(file.path(cfg$out_dir, "metagenomics_vs_mix_scatter.pdf"),
         plot = p_meta, width = 8, height = 7)

  write.table(comp,
              file.path(cfg$out_dir, "metagenomics_vs_mix.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("Metagenomics scatter saved.")

} else {
  msg("Metagenomics file not found (%s) — skipping.",
      cfg$metagenomics_file %||% "NULL")
}


msg("Done. Outputs in: %s/", cfg$out_dir)
