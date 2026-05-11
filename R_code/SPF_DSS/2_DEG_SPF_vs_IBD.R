library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

# Ensure dplyr verbs are not masked by plyr if loaded in the same session
filter    <- dplyr::filter
select    <- dplyr::select
summarise <- dplyr::summarise
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange

dir <- "D:/microbiome_revision/compare_SPF_IBD"
setwd(dir)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
SPF_raw <- read.table("GSM8921752_SPF_stdata.tsv",
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)
IBD_raw <- read.table("GSM8921753_SPFIBD_stdata.tsv",
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)

cat("SPF spots:", nrow(SPF_raw), "| IBD spots:", nrow(IBD_raw), "\n")

# ── 2. Bacterial gene prefixes & species map ─────────────────────────────────
species_map <- c(
  "BcFMB" = "B. choerinum",
  "DY261" = "B. coccoides",
  "BO225" = "D. newyorkensis",
  "EZ315" = "D. freteri",
  "C5O23" = "D. muris",
  "EUAN"  = "E. angustum",
  "KDJ91" = "E. hallii",
  "ARA08" = "E. ramulus",
  "GXM20" = "F. plautii",
  "E0E29" = "I. timonensis",
  "F8244" = "L. gasseri",
  "OCV80" = "O. acetigenes",
  "CRIB"  = "R. ilealis",
  "RGna"  = "R. gnavus",
  "BUB60" = "S. termitidis",
  "GN277" = "S. musculi"
)

prefixes <- names(species_map)
pattern  <- paste0("^(", paste(prefixes, collapse = "|"), ")_RS")

bact_spf <- grep(pattern, colnames(SPF_raw), value = TRUE)
bact_ibd <- grep(pattern, colnames(IBD_raw), value = TRUE)
shared   <- intersect(bact_spf, bact_ibd)

cat("Bacterial genes — SPF:", length(bact_spf),
    "| IBD:", length(bact_ibd),
    "| Shared:", length(shared), "\n")

# ── 3. Build Seurat objects (genes × spots; shared bacterial genes only) ──────
# Seurat expects genes as rows, spots as columns → transpose
spf_mat <- t(as.matrix(SPF_raw[, shared]))
ibd_mat <- t(as.matrix(IBD_raw[, shared]))

spf_obj <- CreateSeuratObject(counts = spf_mat, project = "SPF")
ibd_obj <- CreateSeuratObject(counts = ibd_mat, project = "IBD")

spf_obj$condition <- "SPF"
ibd_obj$condition <- "IBD"

combined <- merge(spf_obj, y = ibd_obj, add.cell.ids = c("SPF", "IBD"))
Idents(combined) <- "condition"

# ── 4. Normalize ──────────────────────────────────────────────────────────────
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
combined <- JoinLayers(combined)

# ── 5. DEG: IBD vs SPF (Wilcoxon, all shared bacterial genes) ─────────────────
# Seurat may sanitize gene names — restrict to genes actually in the object
shared <- rownames(combined)
cat("Running FindMarkers (Wilcoxon) for", length(shared), "bacterial genes...\n")

deg <- FindMarkers(
  combined,
  ident.1         = "IBD",
  ident.2         = "SPF",
  features        = shared,
  test.use        = "wilcox",
  logfc.threshold = 0,        # keep all genes
  min.pct         = 0.001,     # expressed in ≥1% of spots in at least one group
  verbose         = FALSE
)

deg$gene    <- rownames(deg)
deg$prefix  <- sub("[-_]RS.*", "", deg$gene)  # handle both _ and - (Seurat sanitizes _ to -)
deg$species <- species_map[deg$prefix]

# Save full results
write.table(deg, "DEG_IBD_vs_SPF_full.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Full DEG results saved: DEG_IBD_vs_SPF_full.tsv\n")

# ── 6. Per-species summary ────────────────────────────────────────────────────
summary_tbl <- deg %>%
  group_by(prefix, species) %>%
  dplyr::summarise(
    n_tested      = n(),
    n_up_IBD      = sum(p_val_adj < 0.2 & avg_log2FC >  0.1, na.rm = TRUE),
    n_down_IBD    = sum(p_val_adj < 0.2 & avg_log2FC < -0.1, na.rm = TRUE),
    n_sig_total   = sum(p_val_adj < 0.2, na.rm = TRUE),
    median_log2FC = round(median(avg_log2FC, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_total))

print(summary_tbl, n = Inf)
write.table(summary_tbl, "DEG_per_species_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Per-species summary saved: DEG_per_species_summary.tsv\n")

# ── 7. Top significant genes per species ─────────────────────────────────────
top_genes <- deg %>%
  filter(p_val_adj < 0.05) %>%
  group_by(prefix, species) %>%
  slice_max(order_by = abs(avg_log2FC), n = 10) %>%
  arrange(species, desc(abs(avg_log2FC)))

write.table(top_genes, "DEG_top10_per_species.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Top 10 DEGs per species saved: DEG_top10_per_species.tsv\n")

# ── 8. Volcano plot per species ───────────────────────────────────────────────
deg$direction <- "NS"
deg$direction[deg$p_val_adj < 0.05 & deg$avg_log2FC >  0.5] <- "Up in IBD"
deg$direction[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.5] <- "Down in IBD"

# Label top 5 genes per species (by adj p-value)
top_labels <- deg %>%
  filter(p_val_adj < 0.05) %>%
  group_by(species) %>%
  slice_min(order_by = p_val_adj, n = 5) %>%
  ungroup()

p_volcano <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                              color = direction)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 2, max.overlaps = 10, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("NS" = "grey70", "Up in IBD" = "firebrick", "Down in IBD" = "steelblue")
  ) +
  facet_wrap(~species, scales = "free_y", ncol = 4) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.4) +
  theme_bw(base_size = 10) +
  labs(
    title  = "DEG: IBD vs SPF — Bacterial species",
    x      = "avg log2FC (IBD vs SPF)",
    y      = "-log10(adj p-value)",
    color  = NULL
  ) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold.italic"))

ggsave("volcano_per_species.pdf", p_volcano, width = 20, height = 18)
ggsave("volcano_per_species.png", p_volcano, width = 20, height = 18, dpi = 150)
cat("Volcano plot saved: volcano_per_species.pdf / .png\n")

# ── 9. Summary bar plot: n significant DEGs per species ──────────────────────
bar_dat <- deg %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up in IBD", "Down in IBD")) %>%
  count(species, direction) %>%
  mutate(n_plot = ifelse(direction == "Down in IBD", -n, n))

p_bar <- ggplot(bar_dat, aes(x = reorder(species, n_plot), y = n_plot, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up in IBD" = "firebrick", "Down in IBD" = "steelblue")) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  theme_bw(base_size = 25) +
  labs(
    title = "Significant DEGs per bacterial species (|log2FC| > 0.5, adj p < 0.05)",
    x = NULL, y = "Number of DEGs", fill = NULL
  ) +
  theme(
    axis.text.y  = element_text(face = "italic"),
    axis.text.x  = element_text(size = 22),
    axis.title.y = element_text(size = 22)
  )

ggsave("barplot_DEG_counts.pdf", p_bar, width = 10, height = 7)
ggsave("barplot_DEG_counts.png", p_bar, width = 10, height = 7, dpi = 150)
cat("Bar plot saved: barplot_DEG_counts.pdf / .png\n")

# ── 10. Dot plot: top 15 up in IBD + top 15 up in SPF ────────────────────────
top_up_IBD <- deg %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 15) %>%
  dplyr::mutate(group = "Up in IBD")

top_up_SPF <- deg %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  dplyr::slice_min(order_by = avg_log2FC, n = 15) %>%
  dplyr::mutate(group = "Up in SPF")

dot_genes <- bind_rows(top_up_IBD, top_up_SPF)
dot_features <- dot_genes$gene  # ordered: IBD-up first, then SPF-up

# Helper: extract DotPlot data, cap pct.exp at 20, and build ggplot manually
build_dotplot <- function(seurat_obj, features, group_by = "condition",
                          n_up = NULL, title = NULL,
                          base_size = 10, show_legend = TRUE) {
  dat <- DotPlot(seurat_obj, features = features,
                 group.by = group_by)$data
  dat$features.plot <- factor(dat$features.plot, levels = features)
  dat$pct.exp <- pmin(dat$pct.exp, 20)   # cap at 20

  p <- ggplot(dat, aes(x = features.plot, y = id,
                       size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    scale_size_continuous(range = c(1, 8), breaks = c(5, 10, 15, 20),
                          limits = c(0, 20), name = "% Expressed") +
    scale_color_gradient(low = "lightgrey", high = "firebrick",
                         name = "Avg Expression") +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x     = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y     = element_text(face = "bold"),
      legend.position = if (show_legend) "right" else "none"
    ) +
    labs(title = title, x = NULL, y = "Condition")

  if (!is.null(n_up) && n_up > 0 && n_up < length(features)) {
    n_dn <- length(features) - n_up
    p <- p +
      geom_vline(xintercept = n_up + 0.5,
                 linetype = "dashed", color = "grey40", linewidth = 0.4) +
      annotate("text", x = n_up / 2 + 0.5,        y = 2.6,
               label = "Up in IBD", size = 3, color = "firebrick",  fontface = "bold") +
      annotate("text", x = n_up + n_dn / 2 + 0.5, y = 2.6,
               label = "Up in SPF", size = 3, color = "steelblue", fontface = "bold")
  }
  return(p)
}

# Build overall top-15 dot plot
p_dot <- build_dotplot(combined, dot_features,
                       n_up   = nrow(top_up_IBD),
                       title  = "Top 15 DEGs: IBD vs SPF")

ggsave("dotplot_top15_DEGs.pdf", p_dot, width = 14, height = 5)
ggsave("dotplot_top15_DEGs.png", p_dot, width = 14, height = 5, dpi = 150)
cat("Dot plot saved: dotplot_top15_DEGs.pdf / .png\n")

# ── 11. Per-species dot plots — publication quality ───────────────────────────
library(patchwork)

pub_theme <- theme_classic(base_size = 9) +
  theme(
    axis.text.x    = element_text(size = 7, angle = 45, hjust = 1, color = "black"),
    axis.text.y    = element_text(size = 9, face = "bold", color = "black"),
    axis.title     = element_blank(),
    plot.title     = element_text(size = 10, face = "bold.italic", hjust = 0.5,
                                  margin = margin(b = 4)),
    legend.text    = element_text(size = 7),
    legend.title   = element_text(size = 7.5, face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin    = margin(6, 8, 4, 6)
  )

make_pub_dotplot <- function(sp_prefix, sp_name, show_legend = FALSE) {
  sp_up_IBD <- deg %>%
    dplyr::filter(prefix == sp_prefix, p_val_adj < 0.2, avg_log2FC > 0) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 10)

  sp_up_SPF <- deg %>%
    dplyr::filter(prefix == sp_prefix, p_val_adj < 0.2, avg_log2FC < 0) %>%
    dplyr::slice_min(order_by = avg_log2FC, n = 10)

  sp_genes <- bind_rows(
    dplyr::mutate(sp_up_IBD, group = "Up in IBD"),
    dplyr::mutate(sp_up_SPF, group = "Up in SPF")
  )

  if (nrow(sp_genes) < 2) {
    message("Skipping ", sp_name, ": fewer than 2 significant DEGs.")
    return(NULL)
  }

  n_up <- nrow(sp_up_IBD)
  n_dn <- nrow(sp_up_SPF)

  p <- build_dotplot(combined, sp_genes$gene,
                     n_up        = n_up,
                     title       = sp_name,
                     base_size   = 9,
                     show_legend = show_legend) +
    pub_theme +
    labs(title = sp_name)

  if (n_up > 0 && n_dn > 0) {
    p <- p +
      annotate("text", x = n_up / 2 + 0.5, y = 2.52,
               label = "IBD\u2191", size = 2.5,
               color = "#B22222", fontface = "bold") +
      annotate("text", x = n_up + n_dn / 2 + 0.5, y = 2.52,
               label = "SPF\u2191", size = 2.5,
               color = "#2166AC", fontface = "bold")
  }

  return(p)
}

# ── Build plot list ───────────────────────────────────────────────────────────
all_plots <- lapply(names(species_map), function(pfx) {
  make_pub_dotplot(pfx, species_map[[pfx]], show_legend = FALSE)
})
names(all_plots) <- names(species_map)
all_plots <- Filter(Negate(is.null), all_plots)

# Extract shared legend from one representative plot
legend_plot <- make_pub_dotplot(names(all_plots)[1],
                                 species_map[[names(all_plots)[1]]],
                                 show_legend = TRUE)
shared_legend <- cowplot::get_legend(
  legend_plot + theme(legend.position = "right",
                      legend.box = "vertical")
)

# ── Arrange 4 panels per page (2 × 2) with shared legend ─────────────────────
plots_per_page <- 4
n_plots  <- length(all_plots)
n_pages  <- ceiling(n_plots / plots_per_page)

pdf("dotplot_per_species_pub.pdf", width = 14, height = 7,
    useDingbats = FALSE)  # useDingbats=FALSE ensures editable text in Illustrator

for (i in seq_len(n_pages)) {
  idx       <- ((i - 1) * plots_per_page + 1):min(i * plots_per_page, n_plots)
  page_plots <- all_plots[idx]

  # Pad incomplete last page with blank spacers
  while (length(page_plots) < plots_per_page) {
    page_plots <- c(page_plots, list(patchwork::plot_spacer()))
  }

  grid_panel <- patchwork::wrap_plots(page_plots, ncol = 2) +
    patchwork::plot_annotation(
      title   = "DEG Dot Plots: IBD vs SPF",
      caption = paste0("Page ", i, " of ", n_pages,
                       " | adj p < 0.2, top 10 per direction"),
      theme   = theme(
        plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 7, color = "grey50", hjust = 1)
      )
    )

  # Attach shared legend to the right
  print(cowplot::plot_grid(grid_panel, shared_legend,
                           ncol = 2, rel_widths = c(1, 0.12)))
}
dev.off()
cat("Publication-quality per-species dot plots saved: dotplot_per_species_pub.pdf\n")

# ── 12. Spatial expression maps for top 15 DEGs ──────────────────────────────
library(tidyr)

# Top 100 DEGs by adjusted p-value (most significant, any direction)
top15_spatial <- deg %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::slice_min(order_by = p_val_adj, n = 100) %>%
  dplyr::arrange(p_val_adj)
spatial_features <- top15_spatial$gene

# ── Parse coordinates from cell barcodes (format: CONDITION_XxY) ─────────────
cell_meta <- data.frame(
  cell      = colnames(combined),
  condition = combined$condition,
  stringsAsFactors = FALSE
)
xy_part          <- sub("^(SPF|IBD)_", "", cell_meta$cell)
cell_meta$x_coord <- as.numeric(sub("x\\d+$", "", xy_part))
cell_meta$y_coord <- as.numeric(sub("^\\d+x",  "", xy_part))

# ── Extract normalised expression for top 100 genes ──────────────────────────
expr_mat <- GetAssayData(combined, layer = "data")[spatial_features, , drop = FALSE]
expr_df  <- as.data.frame(t(as.matrix(expr_mat)))
expr_df$cell <- rownames(expr_df)

# Merge coordinates with expression
spatial_df <- merge(cell_meta, expr_df, by = "cell")

# ── Repair SPF column 49: replace with average of lanes 48 and 50 ────────────
spf_48 <- spatial_df[spatial_df$condition == "SPF" & spatial_df$x_coord == 48, ]
spf_50 <- spatial_df[spatial_df$condition == "SPF" & spatial_df$x_coord == 50, ]

cat("SPF spots at x=48:", nrow(spf_48), "| x=49 (to be replaced):",
    nrow(spatial_df[spatial_df$condition == "SPF" & spatial_df$x_coord == 49, ]),
    "| x=50:", nrow(spf_50), "\n")

# Remove existing x=49 SPF spots (may be corrupt/empty) and replace with interpolation
spatial_df <- spatial_df[!(spatial_df$condition == "SPF" & spatial_df$x_coord == 49), ]

common_y <- intersect(spf_48$y_coord, spf_50$y_coord)
spf_48_m <- spf_48[spf_48$y_coord %in% common_y, ]
spf_50_m <- spf_50[spf_50$y_coord %in% common_y, ]
spf_48_m <- spf_48_m[order(spf_48_m$y_coord), ]
spf_50_m <- spf_50_m[order(spf_50_m$y_coord), ]

interp <- spf_48_m
interp$cell    <- paste0("SPF_49x", interp$y_coord)
interp$x_coord <- 49
interp[, spatial_features] <- (spf_48_m[, spatial_features] + spf_50_m[, spatial_features]) / 2

spatial_df <- rbind(spatial_df, interp)
cat("SPF column 49 repaired by averaging lanes 48 and 50 (",
    nrow(interp), "spots added)\n")

# Reshape to long format
spatial_long <- tidyr::pivot_longer(
  spatial_df,
  cols      = all_of(spatial_features),
  names_to  = "gene",
  values_to = "expression"
)
spatial_long$gene      <- factor(spatial_long$gene, levels = spatial_features)
spatial_long$condition <- factor(spatial_long$condition, levels = c("SPF", "IBD"))

# Add DEG direction label for facet titles
gene_labels <- setNames(
  paste0(top15_spatial$gene, "  (", ifelse(top15_spatial$avg_log2FC > 0,
         "↑IBD", "↑SPF"), ", padj=", formatC(top15_spatial$p_val_adj,
         format = "e", digits = 1), ")"),
  top15_spatial$gene
)
spatial_long$gene_label <- gene_labels[as.character(spatial_long$gene)]
spatial_long$gene_label <- factor(spatial_long$gene_label,
                                   levels = gene_labels[spatial_features])

# ── Shared expression limits for consistent colour scale ─────────────────────
expr_range <- quantile(spatial_long$expression, c(0.02, 0.98), na.rm = TRUE)

# ── Plot: 4 genes per page, SPF | IBD columns ────────────────────────────────
spatial_theme <- theme_void(base_size = 9) +
  theme(
    strip.text        = element_text(size = 7.5, face = "bold", margin = margin(b = 3)),
    legend.position   = "right",
    legend.title      = element_text(size = 7.5, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.margin       = margin(4, 4, 4, 4)
  )

genes_per_page <- 4
n_spatial_pages <- ceiling(length(spatial_features) / genes_per_page)

pdf("spatial_top100_DEGs.pdf", width = 8, height = 10, useDingbats = FALSE)

for (i in seq_len(n_spatial_pages)) {
  idx         <- ((i - 1) * genes_per_page + 1):min(i * genes_per_page,
                                                      length(spatial_features))
  page_genes  <- spatial_features[idx]
  page_labels <- gene_labels[page_genes]
  page_data   <- spatial_long[spatial_long$gene %in% page_genes, ]
  page_data$gene_label <- factor(page_data$gene_label,
                                  levels = page_labels[page_genes])

  p_spatial <- ggplot(page_data,
                      aes(x = x_coord, y = -y_coord,   # flip Y: origin top-left
                          color = expression)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_gradientn(
      colours = c("grey90", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
      limits  = expr_range,
      oob     = scales::squish,
      name    = "Normalised\nExpression"
    ) +
    facet_grid(gene_label ~ condition, switch = "y") +
    coord_fixed() +
    spatial_theme +
    labs(title = paste0("Spatial Expression — Top DEGs (page ", i,
                        " of ", n_spatial_pages, ")"))

  print(p_spatial)
}

dev.off()
cat("Spatial expression maps saved: spatial_top100_DEGs.pdf\n")

# ── 13. Custom dotplot for selected genes ─────────────────────────────────────
custom_genes <- c(
  "KDJ91-RS14580", "KDJ91-RS05340", "KDJ91-RS07735", "KDJ91-RS10200",
  "RGna-RS11260",  "BUB60-RS19410", "GXM20-RS16595", "GN277-RS23215",
  "GXM20-RS16635", "E0E29-RS05065", "ARA08-RS12880", "KDJ91-RS07740",
  "F8244-RS09290", "KDJ91-RS05345", "GN277-RS26490"
)

# Keep only genes present in the Seurat object
custom_genes <- custom_genes[custom_genes %in% rownames(combined)]
cat("Custom genes found in object:", length(custom_genes), "/", 15, "\n")

p_custom <- build_dotplot(combined, custom_genes,
                          title = "Selected DEGs: IBD vs SPF",
                          base_size = 25) +
  theme(axis.text.x = element_text(size = 22))

ggsave("dotplot_custom_genes.pdf", p_custom, width = 14, height = 4)
ggsave("dotplot_custom_genes.png", p_custom, width = 14, height = 4, dpi = 150)
cat("Custom gene dot plot saved: dotplot_custom_genes.pdf / .png\n")

# ── 14. Spatial plots for custom genes ───────────────────────────────────────
# Filter to genes present in spatial_df (which already has top-100 columns)
custom_spatial_genes <- custom_genes[custom_genes %in% colnames(spatial_df)]

if (length(custom_spatial_genes) == 0) {
  # Need to extract expression for custom genes not yet in spatial_df
  custom_expr_mat <- GetAssayData(combined, layer = "data")[
    custom_genes[custom_genes %in% rownames(combined)], , drop = FALSE]
  custom_expr_df  <- as.data.frame(t(as.matrix(custom_expr_mat)))
  custom_expr_df$cell <- rownames(custom_expr_df)
  custom_spatial_df <- merge(
    spatial_df[, c("cell", "condition", "x_coord", "y_coord")],
    custom_expr_df, by = "cell")
  custom_spatial_genes <- custom_genes[custom_genes %in% rownames(combined)]
} else {
  custom_spatial_df <- spatial_df
}

custom_spatial_long <- tidyr::pivot_longer(
  custom_spatial_df,
  cols      = all_of(custom_spatial_genes),
  names_to  = "gene",
  values_to = "expression"
)
custom_spatial_long$gene      <- factor(custom_spatial_long$gene, levels = custom_spatial_genes)
custom_spatial_long$condition <- factor(custom_spatial_long$condition, levels = c("SPF", "IBD"))

custom_expr_range <- quantile(custom_spatial_long$expression, c(0.02, 0.98), na.rm = TRUE)

# Build one mini-plot per gene (SPF | IBD side by side), then arrange 4 per row
gene_plots <- lapply(custom_spatial_genes, function(g) {
  gdata <- custom_spatial_long[custom_spatial_long$gene == g, ]
  ggplot(gdata, aes(x = x_coord, y = -y_coord, color = expression)) +
    geom_point(size = 0.4, alpha = 0.85) +
    scale_color_gradientn(
      colours = c("grey90", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
      limits  = custom_expr_range,
      oob     = scales::squish,
      name    = "Expr"
    ) +
    facet_wrap(~ condition, nrow = 1) +
    coord_fixed() +
    spatial_theme +
    labs(title = g) +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5),
      strip.text      = element_text(size = 6.5, face = "bold"),
      legend.position = "none",
      plot.margin     = margin(1, 1, 1, 1)
    )
})

# Add shared legend from last gene
legend_gene <- custom_spatial_long[custom_spatial_long$gene == custom_spatial_genes[1], ]
p_legend <- ggplot(legend_gene, aes(x = x_coord, y = -y_coord, color = expression)) +
  geom_point(size = 0.4) +
  scale_color_gradientn(
    colours = c("grey90", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
    limits  = custom_expr_range, oob = scales::squish,
    name    = "Normalised\nExpression"
  ) +
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text  = element_text(size = 6.5),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width  = unit(0.3, "cm"))
shared_legend <- cowplot::get_legend(p_legend)

# Arrange: 4 genes per row (= 8 condition panels per row), all on one page
n_cols  <- 4   # genes per row
n_rows  <- ceiling(length(gene_plots) / n_cols)
grid    <- patchwork::wrap_plots(gene_plots, ncol = n_cols) +
  patchwork::plot_layout(heights = rep(1, n_rows)) &
  theme(plot.margin = margin(1, 1, 1, 1))

# Page size: each gene ~3.5 cm wide × 2 panels, plus legend column
pdf("spatial_custom_genes.pdf",
    width  = n_cols * 3.8 + 1,   # gene columns + legend
    height = n_rows * 2.8,
    useDingbats = FALSE)

print(cowplot::plot_grid(grid, shared_legend,
                         ncol = 2, rel_widths = c(1, 0.06)))
dev.off()
cat("Custom gene spatial plots saved: spatial_custom_genes.pdf\n")

cat("\nDone.\n")
