library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(pheatmap)
library(ggrepel)
library(scales)

dir <- "D:/microbiome_revision/compare_SPF_IBD"
setwd(dir)

# ═════════════════════════════════════════════════════════════════════════════
# SHARED SETUP  (identical to script 4 — kept here for standalone execution)
# ═════════════════════════════════════════════════════════════════════════════

cat("Loading data...\n")
SPF_raw <- read.table("GSM8921752_SPF_stdata.tsv",
                      sep="\t", header=TRUE, row.names=1,
                      stringsAsFactors=FALSE, check.names=FALSE)
IBD_raw <- read.table("GSM8921753_SPFIBD_stdata.tsv",
                      sep="\t", header=TRUE, row.names=1,
                      stringsAsFactors=FALSE, check.names=FALSE)

species_map <- c(
  "BcFMB"="B. choerinum",  "DY261"="B. coccoides",
  "BO225"="D. newyorkensis","EZ315"="D. freteri",
  "C5O23"="D. muris",      "EUAN" ="E. angustum",
  "KDJ91"="E. hallii",     "ARA08"="E. ramulus",
  "GXM20"="F. plautii",    "E0E29"="I. timonensis",
  "F8244"="L. gasseri",    "OCV80"="O. acetigenes",
  "CRIB" ="R. ilealis",    "RGna" ="R. gnavus",
  "BUB60"="S. termitidis", "GN277"="S. musculi"
)
sp_names    <- as.character(species_map)
bact_pat    <- paste0("^(", paste(names(species_map), collapse="|"), ")_RS")
bact_spf    <- grep(bact_pat, colnames(SPF_raw), value=TRUE)
bact_ibd    <- grep(bact_pat, colnames(IBD_raw), value=TRUE)
host_cols   <- intersect(setdiff(colnames(SPF_raw), bact_spf),
                         setdiff(colnames(IBD_raw), bact_ibd))

parse_pos <- function(f) {
  tok <- trimws(unlist(strsplit(paste(readLines(f), collapse=","), ",")))
  tok[grepl("^\\d+x\\d+$", tok)]
}
spf_pos <- parse_pos("SPF_position.txt")
ibd_pos <- parse_pos("IBD_position.txt")

make_spot_df <- function(raw, cond, bcols) {
  m <- data.frame(spot=rownames(raw), condition=cond,
                  x_coord=as.numeric(sub("x\\d+$","",rownames(raw))),
                  y_coord=as.numeric(sub("^\\d+x", "",rownames(raw))),
                  stringsAsFactors=FALSE)
  m$pos_key <- paste0(m$x_coord,"x",m$y_coord)
  for (pfx in names(species_map)) {
    g <- grep(paste0("^",pfx,"_RS"), bcols, value=TRUE)
    m[[species_map[[pfx]]]] <- if(length(g)) rowSums(raw[,g,drop=FALSE],na.rm=TRUE) else 0
  }
  cbind(m, raw[, host_cols, drop=FALSE])
}

spf_df <- make_spot_df(SPF_raw, "SPF", bact_spf)
ibd_df <- make_spot_df(IBD_raw, "IBD", bact_ibd)

# Repair SPF column 49
dcols  <- c(sp_names, host_cols)
s48    <- spf_df[spf_df$x_coord==48,]; s50 <- spf_df[spf_df$x_coord==50,]
spf_df <- spf_df[spf_df$x_coord!=49,]
cy     <- intersect(s48$y_coord, s50$y_coord)
a48    <- s48[match(cy, s48$y_coord),]; a50 <- s50[match(cy, s50$y_coord),]
itp    <- a48
itp$spot <- paste0("49x",cy); itp$x_coord <- 49; itp$pos_key <- paste0("49x",cy)
itp[,dcols] <- (a48[,dcols] + a50[,dcols]) / 2
spf_df <- rbind(spf_df, itp)

# Position filter & merge
spf_df <- spf_df[spf_df$pos_key %in% spf_pos,]
ibd_df <- ibd_df[ibd_df$pos_key %in% ibd_pos,]
all_df <- rbind(spf_df, ibd_df)
cat("Spots — SPF:", nrow(spf_df), "| IBD:", nrow(ibd_df), "\n")

# Normalise
host_raw  <- as.matrix(all_df[, host_cols])
host_norm <- log1p(sweep(host_raw, 1, rowSums(host_raw,na.rm=TRUE)+1, "/")*10000)
bact_mat  <- log1p(as.matrix(all_df[, sp_names]))

# Host modules (case-insensitive match)
mdef <- list(
  Inflammation   =c("Il1b","Tnf","Il6","Cxcl1","Cxcl2","Nos2",
                    "S100a8","S100a9","Ccl2","Il18","Nlrp3","Ptges"),
  Barrier        =c("Muc2","Cldn3","Cldn4","Cldn7","Ocln","Tjp1","Cdh1"),
  AntInflammatory=c("Il10","Tgfb1","Il22","Areg","Il33"),
  Treg           =c("Foxp3","Ctla4","Ikzf2","Il2ra","Tnfrsf18"),
  Goblet         =c("Muc2","Spdef","Atoh1","Klf4","Tff3","Fcgbp"),
  Paneth         =c("Lyz1","Defa1","Defa3","Mmp7","Ang4","Defa17"),
  Enterocyte     =c("Alpi","Fabp1","Fabp2","Slc2a5","Sis","Apoa1"),
  Immune         =c("Ptprc","Cd3e","Cd4","Cd8a","Adgre1","Itgam","Ncam1"),
  Stress         =c("Hspa1a","Atf3","Ddit3","Sqstm1","Hmox1"),
  Proliferation  =c("Mki67","Top2a","Pcna","Mcm2","Mcm6")
)
match_g  <- function(cands, avail) avail[tolower(avail) %in% tolower(cands)]
modules  <- Filter(function(x) length(x)>=2,
                   lapply(mdef, match_g, avail=host_cols))
if (length(modules)==0) stop("No module genes matched. Check host gene name format.")

mod_scores <- as.data.frame(sapply(names(modules), function(m)
  rowMeans(host_norm[, modules[[m]], drop=FALSE], na.rm=TRUE)))

# Fast Spearman via rank + t-approximation
spear <- function(xm, ym) {
  n  <- nrow(xm)
  r  <- cor(apply(xm,2,rank), apply(ym,2,rank), use="pairwise.complete.obs")
  t  <- r * sqrt((n-2)/pmax(1-r^2, 1e-10))
  list(r=r, p=2*pt(-abs(t), df=n-2))
}
idx_s   <- all_df$condition=="SPF"
idx_i   <- all_df$condition=="IBD"
cor_spf <- spear(bact_mat[idx_s,], as.matrix(mod_scores[idx_s,]))
cor_ibd <- spear(bact_mat[idx_i,], as.matrix(mod_scores[idx_i,]))
delta_r <- cor_ibd$r - cor_spf$r
n_tests <- length(cor_spf$r)

# Top significant hits (Bonferroni)
melt_cor <- function(r,p,cond) {
  df <- expand.grid(species=rownames(r), module=colnames(r), stringsAsFactors=FALSE)
  df$r <- as.vector(r); df$p_adj <- as.vector(p)*n_tests; df$cond <- cond; df
}
all_cor  <- rbind(melt_cor(cor_spf$r,cor_spf$p,"SPF"),
                  melt_cor(cor_ibd$r,cor_ibd$p,"IBD"))
top_hits <- all_cor %>% filter(p_adj<0.05) %>%
  arrange(desc(abs(r))) %>% group_by(cond) %>% slice_head(n=6) %>% ungroup()

cat("Setup complete. Building figures...\n\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Bipartite network: species ↔ host modules
#   Edges drawn only for Bonferroni-significant pairs
#   Width ∝ |r|, red = positive, blue = negative
# ═════════════════════════════════════════════════════════════════════════════
make_network_page <- function(r_mat, p_mat, condition) {
  edges <- expand.grid(species=rownames(r_mat), module=colnames(r_mat),
                       stringsAsFactors=FALSE)
  edges$r     <- as.vector(r_mat)
  edges$p_adj <- as.vector(p_mat) * n_tests
  edges        <- edges[edges$p_adj < 0.05, ]
  if (nrow(edges)==0) { message("No sig edges for ", condition); return(NULL) }

  n_sp  <- nrow(r_mat); n_md <- ncol(r_mat)
  sp_nodes  <- data.frame(name=rownames(r_mat), x=0,
                           y=seq(1,0,length.out=n_sp), stringsAsFactors=FALSE)
  mod_nodes <- data.frame(name=colnames(r_mat), x=1,
                           y=seq(1,0,length.out=n_md), stringsAsFactors=FALSE)

  edges <- merge(edges, sp_nodes,  by.x="species", by.y="name")
  names(edges)[names(edges)%in%c("x","y")] <- c("x0","y0")
  edges <- merge(edges, mod_nodes, by.x="module",  by.y="name")
  names(edges)[names(edges)%in%c("x","y")] <- c("x1","y1")
  edges$dir <- ifelse(edges$r>0,"Positive","Negative")

  ggplot() +
    geom_segment(data=edges,
                 aes(x=x0,xend=x1,y=y0,yend=y1,
                     color=dir, linewidth=abs(r)), alpha=0.65) +
    scale_color_manual(values=c(Positive="#B2182B",Negative="#2166AC"),
                       name="Direction") +
    scale_linewidth_continuous(range=c(0.3,3), name="|r|") +
    geom_point(data=sp_nodes,  aes(x=x-0.01,y=y), shape=21, fill="#FDB863",
               color="grey30", size=4) +
    geom_point(data=mod_nodes, aes(x=x+0.01,y=y), shape=22, fill="#B2ABD2",
               color="grey30", size=4) +
    geom_text(data=sp_nodes,  aes(x=x-0.04,y=y, label=name),
              hjust=1, size=2.8, fontface="bold.italic") +
    geom_text(data=mod_nodes, aes(x=x+0.04,y=y, label=name),
              hjust=0, size=2.8, fontface="bold") +
    xlim(-0.7,1.7) +
    theme_void(base_size=10) +
    labs(title=paste0(condition,
                      ": Bacterial–host module association network\n",
                      "(edges = Bonferroni p < 0.05)")) +
    theme(plot.title=element_text(size=10,face="bold",hjust=0.5),
          legend.position="bottom",
          legend.box="horizontal")
}

pdf("host_bacteria_network.pdf", width=11, height=9, useDingbats=FALSE)
pn_spf <- make_network_page(cor_spf$r, cor_spf$p, "SPF")
pn_ibd <- make_network_page(cor_ibd$r, cor_ibd$p, "IBD")
if (!is.null(pn_spf)) print(pn_spf)
if (!is.null(pn_ibd)) print(pn_ibd)
dev.off()
cat("Saved: host_bacteria_network.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — SPF vs IBD correlation scatter
#   Each point = one species × module pair
#   Diagonal = no disease-driven change; deviations = rewiring
# ═════════════════════════════════════════════════════════════════════════════
scat_df <- data.frame(
  species = rep(rownames(cor_spf$r), ncol(cor_spf$r)),
  module  = rep(colnames(cor_spf$r), each=nrow(cor_spf$r)),
  r_spf   = as.vector(cor_spf$r),
  r_ibd   = as.vector(cor_ibd$r),
  p_spf   = as.vector(cor_spf$p)*n_tests,
  p_ibd   = as.vector(cor_ibd$p)*n_tests
)
scat_df$sig   <- scat_df$p_spf<0.05 | scat_df$p_ibd<0.05
scat_df$delta <- scat_df$r_ibd - scat_df$r_spf
scat_df$label <- paste0(scat_df$species,"\n×\n",scat_df$module)

label_thresh <- quantile(abs(scat_df$delta), 0.90, na.rm=TRUE)

p_scat <- ggplot(scat_df, aes(r_spf, r_ibd)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
  geom_hline(yintercept=0, color="grey75", linewidth=0.3) +
  geom_vline(xintercept=0, color="grey75", linewidth=0.3) +
  geom_point(aes(color=delta, size=sig), alpha=0.75) +
  scale_color_gradientn(
    colours=c("#2166AC","white","#B2182B"), limits=c(-0.4,0.4),
    oob=squish, name="Δr\n(IBD−SPF)") +
  scale_size_manual(values=c("FALSE"=1.2,"TRUE"=2.5),
                    labels=c("NS","Sig. in ≥1 condition"), name=NULL) +
  geom_text_repel(
    data=scat_df[abs(scat_df$delta)>=label_thresh & scat_df$sig, ],
    aes(label=label), size=2.1, max.overlaps=20, lineheight=0.85,
    segment.color="grey60", segment.size=0.3) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
           label="Above diagonal = strengthened in IBD\nBelow diagonal = weakened in IBD",
           size=2.8, color="grey40", lineheight=0.9) +
  theme_bw(base_size=10) +
  labs(x="Spearman r  (SPF)", y="Spearman r  (IBD)",
       title="Disease-driven rewiring of host–bacteria associations") +
  theme(legend.position="right",
        plot.title=element_text(size=11,face="bold",hjust=0.5))

ggsave("host_bacteria_cor_scatter.pdf", p_scat,
       width=8, height=7, useDingbats=FALSE)
cat("Saved: host_bacteria_cor_scatter.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Host module spatial expression maps  (SPF | IBD)
#   Each page = one module; consistent colour scale across conditions
# ═════════════════════════════════════════════════════════════════════════════
mod_spatial <- cbind(all_df[,c("condition","x_coord","y_coord")], mod_scores)
mod_long    <- pivot_longer(mod_spatial, cols=names(modules),
                            names_to="module", values_to="score")
mod_long$condition <- factor(mod_long$condition, levels=c("SPF","IBD"))

mod_theme <- theme_void(base_size=9) +
  theme(strip.text=element_text(size=9,face="bold"),
        plot.title=element_text(size=10,face="bold",hjust=0.5),
        legend.key.height=unit(0.5,"cm"), legend.key.width=unit(0.25,"cm"),
        legend.text=element_text(size=6.5), legend.title=element_text(size=7,face="bold"))

pdf("host_module_spatial_maps.pdf", width=9, height=5, useDingbats=FALSE)
for (md in names(modules)) {
  sub  <- mod_long[mod_long$module==md,]
  qlim <- quantile(sub$score, c(0.02,0.98), na.rm=TRUE)
  p <- ggplot(sub, aes(x_coord, -y_coord, color=score)) +
    geom_point(size=0.5, alpha=0.85) +
    scale_color_gradientn(
      colours=c("grey92","#C6DBEF","#4292C6","#08306B"),
      limits=qlim, oob=squish, name="Score") +
    facet_wrap(~condition, nrow=1) +
    coord_fixed() + mod_theme +
    labs(title=paste0(md, " module — spatial expression"))
  print(p)
}
dev.off()
cat("Saved: host_module_spatial_maps.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Spot-level scatter: bacteria log-UMI × module score
#   Coloured by tissue depth (y_coord) to reveal spatial gradient of coupling
# ═════════════════════════════════════════════════════════════════════════════
pdf("host_bacteria_spot_scatter.pdf", width=14, height=5, useDingbats=FALSE)

n_show <- min(nrow(top_hits), 8)
sp_list <- lapply(seq_len(n_show), function(i) {
  h   <- top_hits[i,]
  idx <- all_df$condition == h$cond
  data.frame(
    b     = bact_mat[idx, h$species],
    m     = mod_scores[[h$module]][idx],
    depth = all_df$y_coord[idx],
    xlab  = paste0("log(", h$species, " UMI+1)"),
    ylab  = paste0(h$module, " score"),
    title = paste0("[",h$cond,"]  r=",round(h$r,2))
  )
})

pp_list <- lapply(seq_len(n_show), function(i) {
  d <- sp_list[[i]]
  ggplot(d, aes(b, m, color=depth)) +
    geom_point(size=0.45, alpha=0.6) +
    scale_color_viridis_c(name="Depth\n(y)", direction=-1, option="magma") +
    geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.7) +
    theme_bw(base_size=8) +
    labs(title=d$title[1], x=d$xlab[1], y=d$ylab[1]) +
    theme(plot.title=element_text(size=8,face="bold",hjust=0.5),
          legend.key.height=unit(0.35,"cm"), legend.key.width=unit(0.2,"cm"),
          legend.text=element_text(size=5.5), legend.title=element_text(size=6))
})
print(wrap_plots(pp_list, ncol=4))
dev.off()
cat("Saved: host_bacteria_spot_scatter.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Bacterial Shannon diversity vs host module scores
#   Tests whether spots with richer bacterial communities differ in host state
# ═════════════════════════════════════════════════════════════════════════════
props   <- sweep(bact_mat+1e-9, 1, rowSums(bact_mat+1e-9), "/")
shannon <- -rowSums(props * log(props), na.rm=TRUE)

div_df   <- cbind(data.frame(shannon=shannon, condition=all_df$condition), mod_scores)
div_long <- pivot_longer(div_df, cols=names(modules),
                         names_to="module", values_to="score")
div_long$condition <- factor(div_long$condition, levels=c("SPF","IBD"))

# Per condition × module: add Spearman r annotation
div_ann <- div_long %>%
  group_by(condition, module) %>%
  summarise(
    r   = round(cor(shannon, score, method="spearman", use="complete.obs"), 2),
    p   = cor.test(shannon, score, method="spearman", exact=FALSE)$p.value,
    .groups="drop"
  ) %>%
  mutate(label=paste0("r=",r, ifelse(p<0.05," *","")),
         vjust=ifelse(condition=="SPF", 1.8, 3.2))

p_div <- ggplot(div_long, aes(shannon, score, color=condition)) +
  geom_point(size=0.3, alpha=0.3) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8) +
  geom_text(data=div_ann, aes(label=label, color=condition, vjust=vjust),
            x=Inf, y=Inf, hjust=1.1,
            size=2.6, fontface="bold", show.legend=FALSE) +
  scale_color_manual(values=c(SPF="#4DAF4A", IBD="#E41A1C"), name="Condition") +
  facet_wrap(~module, scales="free_y", ncol=4) +
  theme_bw(base_size=9) +
  labs(x="Bacterial Shannon diversity (per spot)",
       y="Host module score",
       title="Bacterial community diversity vs host gene programs") +
  theme(strip.text=element_text(size=8,face="bold"),
        legend.position="top",
        plot.title=element_text(size=11,face="bold",hjust=0.5))

ggsave("host_bacteria_diversity.pdf", p_div,
       width=14, height=3.6*ceiling(length(modules)/4), useDingbats=FALSE)
cat("Saved: host_bacteria_diversity.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Volcano plots: host genes enriched near each bacterial species
#   Reads the conditional DEG table produced by script 4 if present,
#   otherwise re-runs the test for top-hit species only
# ═════════════════════════════════════════════════════════════════════════════
deg_file <- "host_bacteria_conditional_DEG.tsv"

if (file.exists(deg_file)) {
  cdeg <- read.table(deg_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cat("Loaded conditional DEG from", deg_file, "\n")
} else {
  cat("Re-running conditional DEG for top species...\n")
  gvar  <- apply(host_norm, 2, var, na.rm=TRUE)
  tgene <- names(sort(gvar, decreasing=TRUE))[1:min(2000, ncol(host_norm))]

  run_cdeg <- function(sp, cond, thr=0) {
    idx <- all_df$condition==cond
    prs <- all_df[[sp]][idx] > thr
    if (sum(prs)<10 || sum(!prs)<10) return(NULL)
    h <- host_norm[idx, tgene, drop=FALSE]
    do.call(rbind, lapply(tgene, function(g) {
      w <- suppressWarnings(wilcox.test(h[prs,g], h[!prs,g], exact=FALSE))
      data.frame(gene=g, species=sp, condition=cond,
                 mean_present=mean(h[prs,g], na.rm=TRUE),
                 mean_absent =mean(h[!prs,g],na.rm=TRUE),
                 log2FC=log2((mean(h[prs,g],na.rm=TRUE)+1e-6)/
                             (mean(h[!prs,g],na.rm=TRUE)+1e-6)),
                 p_val=w$p.value, stringsAsFactors=FALSE)
    })) -> res
    res$p_adj <- p.adjust(res$p_val, "BH"); res[order(res$p_adj),]
  }
  cdeg <- do.call(rbind, unlist(lapply(unique(top_hits$species), function(sp)
    lapply(c("SPF","IBD"), function(co) run_cdeg(sp, co))), recursive=FALSE))
  write.table(cdeg, deg_file, sep="\t", quote=FALSE, row.names=FALSE)
}

cdeg <- cdeg %>%
  mutate(direction=case_when(
    p_adj<0.05 & log2FC> 0.5 ~ "Higher near bacteria",
    p_adj<0.05 & log2FC< -0.5~ "Lower near bacteria",
    TRUE                     ~ "NS"),
    panel=paste0(species,"  [",condition,"]"))

top_lab <- cdeg %>% filter(direction!="NS") %>%
  group_by(panel) %>% slice_min(p_adj, n=6) %>% ungroup()

p_volc <- ggplot(cdeg, aes(log2FC, -log10(p_adj+1e-300), color=direction)) +
  geom_point(alpha=0.45, size=0.55) +
  geom_text_repel(data=top_lab, aes(label=gene),
                  size=2.0, max.overlaps=12, show.legend=FALSE,
                  segment.color="grey60", segment.size=0.3) +
  scale_color_manual(values=c("Higher near bacteria"="#B2182B",
                               "Lower near bacteria" ="#2166AC",
                               "NS"="grey75"), name=NULL) +
  geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", alpha=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.4) +
  facet_wrap(~panel, scales="free", ncol=2) +
  theme_bw(base_size=9) +
  labs(x="log2FC  (bacteria-present / absent spots)",
       y="-log10(adj p-value)",
       title="Host gene expression: bacteria-present vs absent spots") +
  theme(strip.text=element_text(size=7.5,face="bold.italic"),
        legend.position="top",
        plot.title=element_text(size=11,face="bold",hjust=0.5))

n_pan <- length(unique(cdeg$panel))
ggsave("host_bacteria_DEG_volcano.pdf", p_volc,
       width=10, height=5*ceiling(n_pan/2), useDingbats=FALSE)
cat("Saved: host_bacteria_DEG_volcano.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Spot heatmap ordered by Inflammation module score
#   Rows = bacterial species; columns = spots (left=most inflamed)
#   Top annotation tracks: Inflammation + Barrier module scores
# ═════════════════════════════════════════════════════════════════════════════
if ("Inflammation" %in% names(modules)) {
  pdf("host_bacteria_spot_heatmap.pdf", width=12, height=5, useDingbats=FALSE)

  for (cond in c("SPF","IBD")) {
    idx <- which(all_df$condition==cond)
    ord <- order(mod_scores$Inflammation[idx], decreasing=TRUE)

    bh  <- t(bact_mat[idx[ord], ])   # species × spots

    # Downsample to ≤600 spots for rendering
    if (ncol(bh)>600) {
      keep <- round(seq(1, ncol(bh), length.out=600))
      bh   <- bh[, keep, drop=FALSE]
      ord2 <- ord[keep]
    } else { ord2 <- ord }

    ann  <- data.frame(
      Inflammation = mod_scores$Inflammation[idx[ord2]],
      row.names    = seq_len(length(ord2))
    )
    if ("Barrier" %in% names(modules))
      ann$Barrier <- mod_scores$Barrier[idx[ord2]]
    if ("AntInflammatory" %in% names(modules))
      ann$AntInflammatory <- mod_scores$AntInflammatory[idx[ord2]]

    colnames(bh) <- seq_len(ncol(bh))

    ann_col_list <- list(
      Inflammation   =colorRampPalette(c("grey95","#B2182B"))(50),
      Barrier        =colorRampPalette(c("grey95","#2166AC"))(50),
      AntInflammatory=colorRampPalette(c("grey95","#1A9641"))(50)
    )
    ann_col_list <- ann_col_list[names(ann_col_list) %in% colnames(ann)]

    pheatmap(bh,
      cluster_rows=TRUE, cluster_cols=FALSE,
      show_colnames=FALSE,
      annotation_col=ann,
      annotation_colors=ann_col_list,
      color=colorRampPalette(c("grey97","#FDD49E","#FC8D59","#D7301F","#7F0000"))(100),
      main=paste0(cond,
                  ": Bacterial species across spots\n",
                  "(columns ordered by Inflammation score, left = most inflamed)"),
      fontsize=8, fontsize_row=7.5, border_color=NA,
      gaps_col=NULL)
  }
  dev.off()
  cat("Saved: host_bacteria_spot_heatmap.pdf\n")
} else {
  cat("Skipping spot heatmap: Inflammation module not found.\n")
}

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 8 — Spatial cluster map
#   k-means (k=6) on z-scored bacterial log-UMI + host module scores per spot
#   Top panel: cluster identity mapped onto tissue grid (SPF | IBD)
#   Bottom panel: cluster profile heatmap — what defines each cluster?
# ═════════════════════════════════════════════════════════════════════════════
set.seed(42)
feat   <- scale(cbind(bact_mat, as.matrix(mod_scores)))
feat[!is.finite(feat)] <- 0
k      <- 6
km     <- kmeans(feat, centers=k, nstart=25, iter.max=150)
all_df$cluster <- factor(km$cluster)

clust_pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628")

# Spatial map
p_cmap <- ggplot(all_df, aes(x_coord, -y_coord, color=cluster)) +
  geom_point(size=0.5, alpha=0.85) +
  scale_color_manual(values=clust_pal[1:k], name="Cluster") +
  facet_wrap(~condition, nrow=1) +
  coord_fixed() +
  theme_void(base_size=9) +
  labs(title=paste0("Spatial clusters (k=",k,
                     "): combined host + bacterial feature space")) +
  theme(plot.title=element_text(size=10,face="bold",hjust=0.5),
        strip.text=element_text(size=9,face="bold"),
        legend.position="right")

# Cluster profile matrix (mean feature per cluster, z-scored columns)
feat_names <- c(sp_names, names(modules))
cpmat <- do.call(rbind, lapply(1:k, function(ci)
  colMeans(cbind(bact_mat, as.matrix(mod_scores))[km$cluster==ci, , drop=FALSE])))
rownames(cpmat) <- paste0("C",1:k)
colnames(cpmat) <- feat_names
cpmat_z <- scale(cpmat)
cpmat_z[!is.finite(cpmat_z)] <- 0

# Cluster size annotation
cl_ann <- data.frame(
  N_spots   = as.integer(table(km$cluster)),
  Condition = sapply(1:k, function(ci)
    names(which.max(table(all_df$condition[km$cluster==ci])))),
  row.names = paste0("C",1:k)
)

pdf("host_bacteria_cluster_spatial.pdf", width=12, height=10, useDingbats=FALSE)
print(p_cmap)
pheatmap(t(cpmat_z),
  annotation_col = cl_ann,
  color   = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  main    = paste0("Cluster profiles — z-scored mean features (k=",k,")"),
  fontsize=8, fontsize_row=7.5, fontsize_col=9,
  border_color="grey85",
  gaps_row=length(sp_names))   # separator between bacteria and module rows
dev.off()
cat("Saved: host_bacteria_cluster_spatial.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 9 — Tissue-zone analysis
#   Divide each section into 4 zones by y_coord quantile
#   (approximates lumen → crypt axis in intestinal ST)
#   Show: bacterial abundance and host module scores per zone, SPF vs IBD
# ═════════════════════════════════════════════════════════════════════════════
zone_df <- all_df[, c("condition","x_coord","y_coord")]
zone_df$zone <- NA

for (cond in c("SPF","IBD")) {
  idx <- zone_df$condition==cond
  yq  <- quantile(zone_df$y_coord[idx], probs=0:4/4, na.rm=TRUE)
  zone_df$zone[idx] <- cut(zone_df$y_coord[idx],
                            breaks=unique(yq), include.lowest=TRUE,
                            labels=paste0("Zone ",seq_len(length(unique(yq))-1)))
}

zone_bact <- cbind(zone_df, bact_mat, mod_scores)
zone_bact_long <- pivot_longer(zone_bact,
  cols=c(sp_names, names(modules)), names_to="feature", values_to="value")
zone_bact_long$type <- ifelse(zone_bact_long$feature %in% sp_names,
                               "Bacteria","Host module")
zone_bact_long$condition <- factor(zone_bact_long$condition, levels=c("SPF","IBD"))

zone_sum <- zone_bact_long %>%
  group_by(condition, zone, feature, type) %>%
  summarise(mean_val=mean(value,na.rm=TRUE), .groups="drop")

# Plot bacteria zones
p_zone_b <- ggplot(zone_sum[zone_sum$type=="Bacteria",],
                   aes(x=zone, y=mean_val, color=condition, group=condition)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2) +
  scale_color_manual(values=c(SPF="#4DAF4A",IBD="#E41A1C"), name=NULL) +
  facet_wrap(~feature, scales="free_y", ncol=4) +
  theme_bw(base_size=8) +
  labs(x="Tissue zone (Zone 1 = top/lumen-proximal)",
       y="Mean log(UMI+1)",
       title="Bacterial species abundance across tissue depth zones") +
  theme(strip.text=element_text(size=7,face="bold.italic"),
        axis.text.x=element_text(size=6,angle=30,hjust=1),
        legend.position="top",
        plot.title=element_text(size=10,face="bold",hjust=0.5))

# Plot module zones
p_zone_m <- ggplot(zone_sum[zone_sum$type=="Host module",],
                   aes(x=zone, y=mean_val, color=condition, group=condition)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2) +
  scale_color_manual(values=c(SPF="#4DAF4A",IBD="#E41A1C"), name=NULL) +
  facet_wrap(~feature, scales="free_y", ncol=4) +
  theme_bw(base_size=8) +
  labs(x="Tissue zone (Zone 1 = top/lumen-proximal)",
       y="Mean module score",
       title="Host module scores across tissue depth zones") +
  theme(strip.text=element_text(size=7,face="bold"),
        axis.text.x=element_text(size=6,angle=30,hjust=1),
        legend.position="top",
        plot.title=element_text(size=10,face="bold",hjust=0.5))

ggsave("host_bacteria_tissue_zones_bacteria.pdf", p_zone_b,
       width=14, height=4*ceiling(length(sp_names)/4), useDingbats=FALSE)
ggsave("host_bacteria_tissue_zones_modules.pdf", p_zone_m,
       width=14, height=4*ceiling(length(modules)/4), useDingbats=FALSE)
cat("Saved: host_bacteria_tissue_zones_bacteria.pdf\n")
cat("Saved: host_bacteria_tissue_zones_modules.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 10 — Ranked correlation lollipop: all species×module pairs per condition
#   Shows the full association landscape sorted by r value
# ═════════════════════════════════════════════════════════════════════════════
loll_df <- all_cor %>%
  mutate(pair  = paste0(species, " × ", module),
         sig   = p_adj < 0.05,
         color = case_when(p_adj<0.05 & r>0 ~ "Positive*",
                           p_adj<0.05 & r<0 ~ "Negative*",
                           TRUE              ~ "NS"))

p_loll <- ggplot(loll_df, aes(x=reorder(pair,r), y=r, color=color)) +
  geom_hline(yintercept=0, linewidth=0.4, color="grey50") +
  geom_segment(aes(xend=pair, yend=0), linewidth=0.5, alpha=0.7) +
  geom_point(aes(size=sig)) +
  scale_color_manual(values=c("Positive*"="#B2182B","Negative*"="#2166AC","NS"="grey65"),
                     name=NULL) +
  scale_size_manual(values=c("TRUE"=2,"FALSE"=1), guide="none") +
  facet_wrap(~cond, ncol=1, scales="free_x") +
  theme_bw(base_size=8) +
  labs(x=NULL, y="Spearman r",
       title="All bacterial species × host module correlations (ranked)") +
  theme(axis.text.x=element_text(size=5.5, angle=60, hjust=1),
        strip.text=element_text(size=9,face="bold"),
        legend.position="top",
        plot.title=element_text(size=10,face="bold",hjust=0.5))

n_pairs <- length(unique(loll_df$pair))
ggsave("host_bacteria_ranked_correlations.pdf", p_loll,
       width=max(10, n_pairs*0.18), height=9, useDingbats=FALSE)
cat("Saved: host_bacteria_ranked_correlations.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
cat("\n═══ All figures complete ═══\n")
cat("1.  host_bacteria_network.pdf\n")
cat("2.  host_bacteria_cor_scatter.pdf\n")
cat("3.  host_module_spatial_maps.pdf\n")
cat("4.  host_bacteria_spot_scatter.pdf\n")
cat("5.  host_bacteria_diversity.pdf\n")
cat("6.  host_bacteria_DEG_volcano.pdf\n")
cat("7.  host_bacteria_spot_heatmap.pdf\n")
cat("8.  host_bacteria_cluster_spatial.pdf\n")
cat("9a. host_bacteria_tissue_zones_bacteria.pdf\n")
cat("9b. host_bacteria_tissue_zones_modules.pdf\n")
cat("10. host_bacteria_ranked_correlations.pdf\n")
