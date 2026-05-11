library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

dir <- "D:/microbiome_revision/compare_SPF_IBD"
setwd(dir)

# ── 1. Load raw data ──────────────────────────────────────────────────────────
cat("Loading data...\n")
SPF_raw <- read.table("GSM8921752_SPF_stdata.tsv",
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)
IBD_raw <- read.table("GSM8921753_SPFIBD_stdata.tsv",
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)

# ── 2. Build Seurat objects (shared bacterial genes) ──────────────────────────
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
shared   <- intersect(
  grep(pattern, colnames(SPF_raw), value = TRUE),
  grep(pattern, colnames(IBD_raw), value = TRUE)
)

spf_obj <- CreateSeuratObject(counts = t(as.matrix(SPF_raw[, shared])), project = "SPF")
ibd_obj <- CreateSeuratObject(counts = t(as.matrix(IBD_raw[, shared])), project = "IBD")
spf_obj$condition <- "SPF"
ibd_obj$condition <- "IBD"

combined <- merge(spf_obj, y = ibd_obj, add.cell.ids = c("SPF", "IBD"))
Idents(combined) <- "condition"
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
combined <- JoinLayers(combined)

# ── 3. Custom genes of interest ───────────────────────────────────────────────
custom_genes <- c(
  "KDJ91-RS14580", "KDJ91-RS05340", "KDJ91-RS07735", "KDJ91-RS10200",
  "RGna-RS11260",  "BUB60-RS19410", "GXM20-RS16595", "GN277-RS23215",
  "GXM20-RS16635", "E0E29-RS05065", "ARA08-RS12880", "KDJ91-RS07740",
  "F8244-RS09290", "KDJ91-RS05345", "GN277-RS26490"
)
custom_genes <- custom_genes[custom_genes %in% rownames(combined)]
cat("Custom genes found in object:", length(custom_genes), "/ 15\n")

# ── 4. Parse useful-spot position files ───────────────────────────────────────
parse_positions <- function(filepath) {
  txt    <- paste(readLines(filepath), collapse = ",")
  tokens <- trimws(unlist(strsplit(txt, ",")))
  tokens[grepl("^\\d+x\\d+$", tokens)]
}

spf_positions <- parse_positions("SPF_position.txt")
ibd_positions <- parse_positions("IBD_position.txt")
cat("Useful spots — SPF:", length(spf_positions),
    "| IBD:", length(ibd_positions), "\n")

# ── 5. Build coordinate table for ALL spots ───────────────────────────────────
cell_meta <- data.frame(
  cell      = colnames(combined),
  condition = combined$condition,
  stringsAsFactors = FALSE
)
xy_part           <- sub("^(SPF|IBD)_", "", cell_meta$cell)
cell_meta$x_coord <- as.numeric(sub("x\\d+$", "", xy_part))
cell_meta$y_coord <- as.numeric(sub("^\\d+x",  "", xy_part))
cell_meta$pos_key <- paste0(cell_meta$x_coord, "x", cell_meta$y_coord)

# ── 6. Extract expression for custom genes (all spots) ────────────────────────
expr_mat <- GetAssayData(combined, layer = "data")[custom_genes, , drop = FALSE]
expr_df  <- as.data.frame(t(as.matrix(expr_mat)))
expr_df$cell <- rownames(expr_df)

all_spatial_df <- merge(
  cell_meta[, c("cell", "condition", "x_coord", "y_coord", "pos_key")],
  expr_df, by = "cell"
)

# ── 7. Repair SPF column 49: replace with average of columns 48 and 50 ────────
spf_48 <- all_spatial_df[all_spatial_df$condition == "SPF" & all_spatial_df$x_coord == 48, ]
spf_50 <- all_spatial_df[all_spatial_df$condition == "SPF" & all_spatial_df$x_coord == 50, ]

cat("SPF spots — x=48:", nrow(spf_48), "| x=49 (defective, removing):",
    nrow(all_spatial_df[all_spatial_df$condition == "SPF" & all_spatial_df$x_coord == 49, ]),
    "| x=50:", nrow(spf_50), "\n")

# Remove defective x=49 SPF spots
all_spatial_df <- all_spatial_df[
  !(all_spatial_df$condition == "SPF" & all_spatial_df$x_coord == 49), ]

# Interpolate: average x=48 and x=50 for rows with matching y
common_y <- intersect(spf_48$y_coord, spf_50$y_coord)
spf_48_m <- spf_48[spf_48$y_coord %in% common_y, ]
spf_50_m <- spf_50[spf_50$y_coord %in% common_y, ]
spf_48_m <- spf_48_m[order(spf_48_m$y_coord), ]
spf_50_m <- spf_50_m[order(spf_50_m$y_coord), ]

interp            <- spf_48_m
interp$cell       <- paste0("SPF_49x", interp$y_coord)
interp$x_coord    <- 49
interp$pos_key    <- paste0("49x", interp$y_coord)
interp[, custom_genes] <- (spf_48_m[, custom_genes] + spf_50_m[, custom_genes]) / 2

all_spatial_df <- rbind(all_spatial_df, interp)
cat("SPF column 49 repaired:", nrow(interp), "interpolated spots added\n")

# ── 8. Apply position filter ──────────────────────────────────────────────────
custom_spatial_df <- all_spatial_df[
  (all_spatial_df$condition == "SPF" & all_spatial_df$pos_key %in% spf_positions) |
  (all_spatial_df$condition == "IBD" & all_spatial_df$pos_key %in% ibd_positions), ]

cat("Spots retained after filtering — SPF:", sum(custom_spatial_df$condition == "SPF"),
    "| IBD:", sum(custom_spatial_df$condition == "IBD"), "\n")

# ── 9. Pivot to long format ───────────────────────────────────────────────────
custom_spatial_long <- tidyr::pivot_longer(
  custom_spatial_df,
  cols      = all_of(custom_genes),
  names_to  = "gene",
  values_to = "expression"
)
custom_spatial_long$gene      <- factor(custom_spatial_long$gene, levels = custom_genes)
custom_spatial_long$condition <- factor(custom_spatial_long$condition, levels = c("SPF", "IBD"))

custom_expr_range <- quantile(custom_spatial_long$expression, c(0.02, 0.98), na.rm = TRUE)

# ── 10. Build per-gene spatial plots ──────────────────────────────────────────
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

gene_plots <- lapply(custom_genes, function(g) {
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

# ── 11. Shared legend ─────────────────────────────────────────────────────────
p_legend <- ggplot(
  custom_spatial_long[custom_spatial_long$gene == custom_genes[1], ],
  aes(x = x_coord, y = -y_coord, color = expression)
) +
  geom_point(size = 0.4) +
  scale_color_gradientn(
    colours = c("grey90", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
    limits  = custom_expr_range, oob = scales::squish,
    name    = "Normalised\nExpression"
  ) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 7, face = "bold"),
    legend.text       = element_text(size = 6.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.3, "cm")
  )
shared_legend <- cowplot::get_legend(p_legend)

# ── 12. Arrange and save PDF ──────────────────────────────────────────────────
n_cols <- 4
n_rows <- ceiling(length(gene_plots) / n_cols)

grid <- patchwork::wrap_plots(gene_plots, ncol = n_cols) +
  patchwork::plot_layout(heights = rep(1, n_rows)) &
  theme(plot.margin = margin(1, 1, 1, 1))

pdf("spatial_custom_genes - Copy.pdf",
    width  = n_cols * 3.8 + 1,
    height = n_rows * 2.8,
    useDingbats = FALSE)

print(cowplot::plot_grid(grid, shared_legend,
                         ncol = 2, rel_widths = c(1, 0.06)))
dev.off()

cat("Saved: spatial_custom_genes - Copy.pdf\n")
cat("Done.\n")
