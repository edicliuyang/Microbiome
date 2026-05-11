library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(pheatmap)

dir <- "D:/microbiome_revision/compare_SPF_IBD"
setwd(dir)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load raw data
# ─────────────────────────────────────────────────────────────────────────────
cat("Loading data...\n")
SPF_raw <- read.table("GSM8921752_SPF_stdata.tsv",
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)
IBD_raw <- read.table("GSM8921753_SPFIBD_stdata.tsv",
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Separate bacterial vs host columns
# ─────────────────────────────────────────────────────────────────────────────
species_map <- c(
  "BcFMB" = "B. choerinum",  "DY261" = "B. coccoides",
  "BO225" = "D. newyorkensis","EZ315" = "D. freteri",
  "C5O23" = "D. muris",      "EUAN"  = "E. angustum",
  "KDJ91" = "E. hallii",     "ARA08" = "E. ramulus",
  "GXM20" = "F. plautii",    "E0E29" = "I. timonensis",
  "F8244" = "L. gasseri",    "OCV80" = "O. acetigenes",
  "CRIB"  = "R. ilealis",    "RGna"  = "R. gnavus",
  "BUB60" = "S. termitidis", "GN277" = "S. musculi"
)

bact_pattern <- paste0("^(", paste(names(species_map), collapse = "|"), ")_RS")
bact_spf  <- grep(bact_pattern, colnames(SPF_raw), value = TRUE)
bact_ibd  <- grep(bact_pattern, colnames(IBD_raw), value = TRUE)
host_spf  <- setdiff(colnames(SPF_raw), bact_spf)
host_ibd  <- setdiff(colnames(IBD_raw), bact_ibd)
host_cols <- intersect(host_spf, host_ibd)

cat("Bacterial genes — SPF:", length(bact_spf), "| IBD:", length(bact_ibd), "\n")
cat("Host genes (shared):", length(host_cols), "\n")
cat("First 10 host column names:", paste(head(host_cols, 10), collapse = ", "), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Parse useful-spot position files
# ─────────────────────────────────────────────────────────────────────────────
parse_positions <- function(filepath) {
  txt    <- paste(readLines(filepath), collapse = ",")
  tokens <- trimws(unlist(strsplit(txt, ",")))
  tokens[grepl("^\\d+x\\d+$", tokens)]
}
spf_positions <- parse_positions("SPF_position.txt")
ibd_positions <- parse_positions("IBD_position.txt")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Build per-spot data frames (meta + per-species UMI + host counts)
# ─────────────────────────────────────────────────────────────────────────────
sp_names <- as.character(species_map)  # italic species name labels

make_spot_df <- function(raw_df, condition, bact_cols) {
  meta <- data.frame(
    spot      = rownames(raw_df),
    condition = condition,
    x_coord   = as.numeric(sub("x\\d+$", "", rownames(raw_df))),
    y_coord   = as.numeric(sub("^\\d+x",  "", rownames(raw_df))),
    stringsAsFactors = FALSE
  )
  meta$pos_key <- paste0(meta$x_coord, "x", meta$y_coord)

  # Sum UMI across all genes per species
  for (pfx in names(species_map)) {
    sp_g <- grep(paste0("^", pfx, "_RS"), bact_cols, value = TRUE)
    meta[[species_map[[pfx]]]] <- if (length(sp_g) > 0)
      rowSums(raw_df[, sp_g, drop = FALSE], na.rm = TRUE) else 0
  }

  # Attach host counts (numeric only)
  host_sub <- raw_df[, host_cols, drop = FALSE]
  cbind(meta, host_sub)
}

spf_df <- make_spot_df(SPF_raw, "SPF", bact_spf)
ibd_df <- make_spot_df(IBD_raw, "IBD", bact_ibd)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Repair SPF column 49 (average of columns 48 and 50)
# ─────────────────────────────────────────────────────────────────────────────
data_cols <- c(sp_names, host_cols)

spf_48 <- spf_df[spf_df$x_coord == 48, ]
spf_50 <- spf_df[spf_df$x_coord == 50, ]
spf_df  <- spf_df[spf_df$x_coord != 49, ]

common_y  <- intersect(spf_48$y_coord, spf_50$y_coord)
s48 <- spf_48[match(common_y, spf_48$y_coord), ]
s50 <- spf_50[match(common_y, spf_50$y_coord), ]

interp              <- s48
interp$spot         <- paste0("49x", common_y)
interp$x_coord      <- 49
interp$pos_key      <- paste0("49x", common_y)
interp[, data_cols] <- (s48[, data_cols] + s50[, data_cols]) / 2

spf_df <- rbind(spf_df, interp)
cat("SPF column 49 repaired:", nrow(interp), "spots added\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Filter to useful positions and merge
# ─────────────────────────────────────────────────────────────────────────────
spf_df <- spf_df[spf_df$pos_key %in% spf_positions, ]
ibd_df <- ibd_df[ibd_df$pos_key %in% ibd_positions, ]
all_df <- rbind(spf_df, ibd_df)
cat("After filtering — SPF:", nrow(spf_df), "| IBD:", nrow(ibd_df), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Normalize host genes (log-CPM) and log1p bacterial UMI
# ─────────────────────────────────────────────────────────────────────────────
host_raw  <- as.matrix(all_df[, host_cols])
lib_size  <- rowSums(host_raw, na.rm = TRUE)
host_norm <- log1p(sweep(host_raw, 1, lib_size + 1, "/") * 10000)

bact_mat  <- log1p(as.matrix(all_df[, sp_names]))

# ─────────────────────────────────────────────────────────────────────────────
# 8. Define host modules and compute scores
#    Module genes matched case-insensitively against available host columns
# ─────────────────────────────────────────────────────────────────────────────
modules_def <- list(
  Inflammation    = c("Il1b","Tnf","Il6","Cxcl1","Cxcl2","Nos2",
                      "S100a8","S100a9","Ccl2","Il18","Nlrp3","Ptges"),
  Barrier         = c("Muc2","Cldn3","Cldn4","Cldn7","Ocln","Tjp1","Cdh1"),
  AntInflammatory = c("Il10","Tgfb1","Il22","Areg","Il33"),
  Treg            = c("Foxp3","Ctla4","Ikzf2","Il2ra","Tnfrsf18"),
  Goblet          = c("Muc2","Spdef","Atoh1","Klf4","Tff3","Fcgbp"),
  Paneth          = c("Lyz1","Defa1","Defa3","Mmp7","Ang4","Defa17"),
  Enterocyte      = c("Alpi","Fabp1","Fabp2","Slc2a5","Sis","Apoa1"),
  Immune          = c("Ptprc","Cd3e","Cd4","Cd8a","Adgre1","Itgam","Ncam1"),
  Stress          = c("Hspa1a","Atf3","Ddit3","Sqstm1","Hmox1"),
  Proliferation   = c("Mki67","Top2a","Pcna","Mcm2","Mcm6")
)

match_genes <- function(candidates, available) {
  available[tolower(available) %in% tolower(candidates)]
}

modules <- lapply(modules_def, match_genes, available = host_cols)
modules <- modules[sapply(modules, length) >= 2]

if (length(modules) == 0) {
  stop(
    "No module genes matched host columns. ",
    "Check whether host gene names are symbols (e.g. Il1b) or ENSEMBL IDs.\n",
    "First 20 host columns: ", paste(head(host_cols, 20), collapse = ", ")
  )
}

cat("\nModules retained (>= 2 genes):\n")
for (m in names(modules))
  cat("  ", m, ":", length(modules[[m]]), "—", paste(modules[[m]], collapse = ", "), "\n")

module_scores <- as.data.frame(sapply(names(modules), function(m)
  rowMeans(host_norm[, modules[[m]], drop = FALSE], na.rm = TRUE)))
rownames(module_scores) <- rownames(all_df)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Spearman correlation: bacterial species × host modules
#    Fast: rank-transform → Pearson; p via t-approximation
# ─────────────────────────────────────────────────────────────────────────────
spear_cor <- function(x_mat, y_mat) {
  xr <- apply(x_mat, 2, rank, na.last = "keep")
  yr <- apply(y_mat, 2, rank, na.last = "keep")
  n  <- nrow(xr)
  r  <- cor(xr, yr, use = "pairwise.complete.obs")
  t  <- r * sqrt((n - 2) / pmax(1 - r^2, 1e-10))
  p  <- 2 * pt(-abs(t), df = n - 2)
  list(r = r, p = p)
}

idx_spf <- all_df$condition == "SPF"
idx_ibd <- all_df$condition == "IBD"

cat("Computing SPF correlations...\n")
cor_spf <- spear_cor(bact_mat[idx_spf, ], as.matrix(module_scores[idx_spf, ]))
cat("Computing IBD correlations...\n")
cor_ibd <- spear_cor(bact_mat[idx_ibd, ], as.matrix(module_scores[idx_ibd, ]))

delta_r <- cor_ibd$r - cor_spf$r   # positive = association strengthened in IBD

n_tests <- length(cor_spf$r)

stars_overlay <- function(p_mat, n_tests) {
  adj <- p_mat * n_tests
  m   <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat),
                dimnames = dimnames(p_mat))
  m[adj < 0.05]  <- "*"
  m[adj < 0.01]  <- "**"
  m[adj < 0.001] <- "***"
  m
}

# ─────────────────────────────────────────────────────────────────────────────
# 10. OUTPUT A: Spearman correlation heatmaps (SPF | IBD | delta)
# ─────────────────────────────────────────────────────────────────────────────
r_lim       <- max(abs(range(c(cor_spf$r, cor_ibd$r), na.rm = TRUE)), 0.1)
cor_breaks  <- seq(-r_lim, r_lim, length.out = 101)
cor_colors  <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

d_lim       <- max(abs(range(delta_r, na.rm = TRUE)), 0.05)
delta_breaks <- seq(-d_lim, d_lim, length.out = 101)
delta_colors <- colorRampPalette(c("#4DAC26", "white", "#D01C8B"))(100)

pdf("host_bacteria_correlation_heatmap.pdf", width = 9, height = 7,
    useDingbats = FALSE)

pheatmap(cor_spf$r,
  color = cor_colors, breaks = cor_breaks,
  main  = "SPF: Bacterial species × Host modules (Spearman r)",
  display_numbers = stars_overlay(cor_spf$p, n_tests),
  number_color = "black", fontsize = 9, border_color = "grey80",
  cluster_rows = TRUE, cluster_cols = TRUE
)

pheatmap(cor_ibd$r,
  color = cor_colors, breaks = cor_breaks,
  main  = "IBD: Bacterial species × Host modules (Spearman r)",
  display_numbers = stars_overlay(cor_ibd$p, n_tests),
  number_color = "black", fontsize = 9, border_color = "grey80",
  cluster_rows = TRUE, cluster_cols = TRUE
)

pheatmap(delta_r,
  color = delta_colors, breaks = delta_breaks,
  main  = "ΔIBD − SPF: shift in bacterial × host association\n(positive = strengthened in IBD)",
  display_numbers = round(delta_r, 2),
  number_color = "black", fontsize = 9, border_color = "grey80",
  cluster_rows = TRUE, cluster_cols = TRUE
)

dev.off()
cat("Saved: host_bacteria_correlation_heatmap.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# 11. Identify top hits (largest |r|, Bonferroni-significant) per condition
# ─────────────────────────────────────────────────────────────────────────────
melt_cor <- function(r_mat, p_mat, cond, n_tests) {
  df <- expand.grid(species = rownames(r_mat), module = colnames(r_mat),
                    stringsAsFactors = FALSE)
  df$r     <- as.vector(r_mat)
  df$p_adj <- as.vector(p_mat) * n_tests
  df$cond  <- cond
  df
}

all_cor <- rbind(
  melt_cor(cor_spf$r, cor_spf$p, "SPF", n_tests),
  melt_cor(cor_ibd$r, cor_ibd$p, "IBD", n_tests)
)

top_hits <- all_cor %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::arrange(desc(abs(r))) %>%
  dplyr::group_by(cond) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::ungroup()

cat("\nTop host-bacteria correlations:\n")
print(as.data.frame(top_hits[, c("cond","species","module","r","p_adj")]))

# ─────────────────────────────────────────────────────────────────────────────
# 12. OUTPUT B: Bivariate spatial maps
#     Three panels per hit: bacteria alone | host module alone | bivariate overlay
#     Bivariate colour: red = bacteria, blue = host, purple = co-high
# ─────────────────────────────────────────────────────────────────────────────
make_bivar_page <- function(sp_name, mod_name, condition, r_val) {
  idx   <- all_df$condition == condition
  b_raw <- log1p(all_df[[sp_name]][idx])
  m_raw <- module_scores[[mod_name]][idx]
  x_c   <- all_df$x_coord[idx]
  y_c   <- -all_df$y_coord[idx]

  scale01 <- function(v) (v - min(v, na.rm=TRUE)) /
                          (max(v, na.rm=TRUE) - min(v, na.rm=TRUE) + 1e-9)
  bivar_col <- rgb(scale01(b_raw), 0, scale01(m_raw))

  base_df <- data.frame(x = x_c, y = y_c, b = b_raw, m = m_raw,
                         bivar = bivar_col)

  grad_bact <- scale_color_gradientn(
    colours = c("grey92","#FDD49E","#FC8D59","#D7301F","#7F0000"),
    name    = "log(UMI+1)", limits = range(b_raw, na.rm = TRUE))

  grad_host <- scale_color_gradientn(
    colours = c("grey92","#C6DBEF","#4292C6","#08306B"),
    name    = "Module score", limits = range(m_raw, na.rm = TRUE))

  pt_size <- 0.5

  p1 <- ggplot(base_df, aes(x, y, color = b)) +
    geom_point(size = pt_size, alpha = 0.9) + grad_bact + coord_fixed() +
    theme_void(base_size = 8) +
    labs(title = paste0("Bacteria: ", sp_name)) +
    theme(plot.title = element_text(size = 7.5, face = "bold.italic", hjust = 0.5),
          legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.2,"cm"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 6))

  p2 <- ggplot(base_df, aes(x, y, color = m)) +
    geom_point(size = pt_size, alpha = 0.9) + grad_host + coord_fixed() +
    theme_void(base_size = 8) +
    labs(title = paste0("Host: ", mod_name)) +
    theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
          legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.2,"cm"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 6))

  p3 <- ggplot(base_df, aes(x, y)) +
    geom_point(color = base_df$bivar, size = pt_size, alpha = 0.9) +
    coord_fixed() + theme_void(base_size = 8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = "red=bacteria | blue=host | purple=co-high",
             size = 2.2, color = "grey30") +
    labs(title = "Bivariate overlay") +
    theme(plot.title = element_text(size = 7.5, hjust = 0.5))

  patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
    patchwork::plot_annotation(
      title = paste0("[", condition, "]  ", sp_name,
                     "  ×  ", mod_name,
                     "  (Spearman r = ", round(r_val, 3), ")"),
      theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
}

pdf("host_bacteria_bivariate_spatial.pdf", width = 13, height = 4.5,
    useDingbats = FALSE)

for (i in seq_len(nrow(top_hits))) {
  h <- top_hits[i, ]
  tryCatch(
    print(make_bivar_page(h$species, h$module, h$cond, h$r)),
    error = function(e)
      message("Skipping ", h$species, " × ", h$module, ": ", e$message)
  )
}

dev.off()
cat("Saved: host_bacteria_bivariate_spatial.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# 13. OUTPUT C: Module score by bacterial abundance quartile
#     One panel per top hit; SPF and IBD shown side by side within each panel
# ─────────────────────────────────────────────────────────────────────────────
quart_rows <- list()

for (i in seq_len(nrow(top_hits))) {
  h   <- top_hits[i, ]
  idx <- all_df$condition == h$cond
  b   <- log1p(all_df[[h$species]][idx])
  m   <- module_scores[[h$module]][idx]

  qs <- quantile(b, probs = 0:4/4, na.rm = TRUE)
  # Force unique breaks if many zeros
  if (length(unique(qs)) < 5) qs <- c(qs[1], quantile(b[b > 0], probs = 1:3/3, na.rm = TRUE), qs[5])
  q_label <- cut(b, breaks = unique(qs), include.lowest = TRUE,
                 labels = paste0("Q", seq_len(length(unique(qs)) - 1)))

  quart_rows[[i]] <- data.frame(
    quartile = q_label,
    score    = m,
    label    = paste0(h$species, " × ", h$module),
    cond     = h$cond,
    stringsAsFactors = FALSE
  )
}

quart_df <- do.call(rbind, quart_rows)
quart_df <- quart_df[!is.na(quart_df$quartile), ]

p_quart <- ggplot(quart_df, aes(x = quartile, y = score, fill = cond)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.4, width = 0.6) +
  scale_fill_manual(values = c("SPF" = "#4DAF4A", "IBD" = "#E41A1C"),
                    name = "Condition") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 9) +
  labs(
    x     = "Bacterial abundance quartile (log UMI+1)",
    y     = "Host module score",
    title = "Host module score across bacterial abundance quartiles"
  ) +
  theme(
    strip.text      = element_text(size = 8, face = "bold"),
    legend.position = "top",
    plot.title      = element_text(size = 10, face = "bold", hjust = 0.5)
  )

n_panels <- length(unique(quart_df$label))
ggsave("host_bacteria_quartile_boxplots.pdf", p_quart,
       width = 12, height = 3.8 * ceiling(n_panels / 3),
       useDingbats = FALSE)
cat("Saved: host_bacteria_quartile_boxplots.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# 14. OUTPUT D: Conditional host DEG — bacteria-present vs absent spots
#     Run for each species appearing in top_hits, per condition
#     Tests top 2000 most-variable host genes (Wilcoxon, BH-adjusted)
# ─────────────────────────────────────────────────────────────────────────────
gene_var  <- apply(host_norm, 2, var, na.rm = TRUE)
test_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(2000, ncol(host_norm))]

run_cond_deg <- function(sp_name, condition, threshold = 0) {
  idx     <- all_df$condition == condition
  present <- all_df[[sp_name]][idx] > threshold

  if (sum(present) < 10 || sum(!present) < 10) {
    message("  Skipping ", sp_name, " [", condition, "]: too few spots in one group")
    return(NULL)
  }

  cat("  Conditional DEG:", sp_name, "[", condition, "] —",
      sum(present), "present vs", sum(!present), "absent spots\n")

  h_sub <- host_norm[idx, test_genes, drop = FALSE]

  res <- do.call(rbind, lapply(test_genes, function(g) {
    x1 <- h_sub[present,  g]
    x0 <- h_sub[!present, g]
    w  <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE))
    data.frame(
      gene         = g,
      species      = sp_name,
      condition    = condition,
      n_present    = sum(present),
      n_absent     = sum(!present),
      mean_present = mean(x1, na.rm = TRUE),
      mean_absent  = mean(x0, na.rm = TRUE),
      log2FC       = log2((mean(x1, na.rm=TRUE) + 1e-6) /
                          (mean(x0, na.rm=TRUE) + 1e-6)),
      p_val        = w$p.value,
      stringsAsFactors = FALSE
    )
  }))

  res$p_adj <- p.adjust(res$p_val, method = "BH")
  res[order(res$p_adj), ]
}

species_to_test <- unique(top_hits$species)
cond_deg_list   <- list()

for (sp in species_to_test) {
  for (cond in c("SPF", "IBD")) {
    res <- run_cond_deg(sp, cond)
    if (!is.null(res)) cond_deg_list[[paste(sp, cond)]] <- res
  }
}

all_cond_deg <- do.call(rbind, cond_deg_list)
write.table(all_cond_deg, "host_bacteria_conditional_DEG.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Conditional DEG table saved: host_bacteria_conditional_DEG.tsv\n")

# Plot: top 10 host DEGs per species × condition
deg_sig <- all_cond_deg %>%
  dplyr::filter(p_adj < 0.05, abs(log2FC) > 0.5) %>%
  dplyr::mutate(
    direction = ifelse(log2FC > 0,
                       "Higher near bacteria", "Lower near bacteria"),
    label     = paste0(species, "  [", condition, "]")
  ) %>%
  dplyr::group_by(label) %>%
  dplyr::slice_max(order_by = abs(log2FC), n = 10) %>%
  dplyr::ungroup()

if (nrow(deg_sig) > 0) {
  p_cdeg <- ggplot(deg_sig,
                   aes(x = reorder(gene, log2FC), y = log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Higher near bacteria" = "#B2182B",
                                  "Lower near bacteria"  = "#2166AC")) +
    facet_wrap(~ label, scales = "free", ncol = 2) +
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(
      title = "Host gene expression: bacteria-present vs absent spots",
      x = NULL, y = "log2FC (present / absent)", fill = NULL
    ) +
    theme(
      strip.text      = element_text(size = 8, face = "bold.italic"),
      legend.position = "top",
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5)
    )

  n_panels <- length(unique(deg_sig$label))
  ggsave("host_bacteria_conditional_DEG.pdf", p_cdeg,
         width = 10, height = 4 * ceiling(n_panels / 2),
         useDingbats = FALSE)
  cat("Saved: host_bacteria_conditional_DEG.pdf\n")
} else {
  cat("No significant conditional DEGs (p_adj < 0.05, |log2FC| > 0.5).\n",
      "Consider relaxing thresholds in the filter above.\n")
}

cat("\nAll host-bacteria co-localization analyses complete.\n")
cat("Outputs:\n")
cat("  host_bacteria_correlation_heatmap.pdf  — Spearman r (SPF, IBD, delta)\n")
cat("  host_bacteria_bivariate_spatial.pdf    — spatial co-localization maps\n")
cat("  host_bacteria_quartile_boxplots.pdf    — module scores by bacterial quartile\n")
cat("  host_bacteria_conditional_DEG.pdf      — host DEGs near bacteria\n")
cat("  host_bacteria_conditional_DEG.tsv      — full conditional DEG table\n")
