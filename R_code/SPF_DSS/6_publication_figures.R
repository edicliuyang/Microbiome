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
# SHARED SETUP
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
sp_names  <- as.character(species_map)
bact_pat  <- paste0("^(", paste(names(species_map), collapse="|"), ")_RS")
bact_spf  <- grep(bact_pat, colnames(SPF_raw), value=TRUE)
bact_ibd  <- grep(bact_pat, colnames(IBD_raw), value=TRUE)
host_cols <- intersect(setdiff(colnames(SPF_raw), bact_spf),
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
    m[[species_map[[pfx]]]] <- if(length(g))
      rowSums(raw[,g,drop=FALSE],na.rm=TRUE) else 0
  }
  cbind(m, raw[, host_cols, drop=FALSE])
}

spf_df <- make_spot_df(SPF_raw,"SPF",bact_spf)
ibd_df <- make_spot_df(IBD_raw,"IBD",bact_ibd)

dcols  <- c(sp_names, host_cols)
s48    <- spf_df[spf_df$x_coord==48,]; s50 <- spf_df[spf_df$x_coord==50,]
spf_df <- spf_df[spf_df$x_coord!=49,]
cy     <- intersect(s48$y_coord, s50$y_coord)
a48    <- s48[match(cy,s48$y_coord),]; a50 <- s50[match(cy,s50$y_coord),]
itp    <- a48; itp$spot <- paste0("49x",cy)
itp$x_coord <- 49; itp$pos_key <- paste0("49x",cy)
itp[,dcols] <- (a48[,dcols]+a50[,dcols])/2
spf_df <- rbind(spf_df,itp)

spf_df <- spf_df[spf_df$pos_key %in% spf_pos,]
ibd_df <- ibd_df[ibd_df$pos_key %in% ibd_pos,]
all_df <- rbind(spf_df,ibd_df)
cat("Spots — SPF:", nrow(spf_df), "| IBD:", nrow(ibd_df), "\n")

host_raw  <- as.matrix(all_df[,host_cols])
host_norm <- log1p(sweep(host_raw,1,rowSums(host_raw,na.rm=TRUE)+1,"/")*10000)
bact_mat  <- log1p(as.matrix(all_df[,sp_names]))

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
match_g <- function(cands,avail) avail[tolower(avail)%in%tolower(cands)]
modules <- Filter(function(x) length(x)>=2, lapply(mdef,match_g,avail=host_cols))
if(length(modules)==0) stop("No module genes matched host columns.")

mod_scores <- as.data.frame(sapply(names(modules), function(m)
  rowMeans(host_norm[,modules[[m]],drop=FALSE],na.rm=TRUE)))

spear <- function(xm,ym) {
  n <- nrow(xm)
  r <- cor(apply(xm,2,rank),apply(ym,2,rank),use="pairwise.complete.obs")
  t <- r*sqrt((n-2)/pmax(1-r^2,1e-10))
  list(r=r, p=2*pt(-abs(t),df=n-2))
}
idx_s   <- all_df$condition=="SPF"
idx_i   <- all_df$condition=="IBD"
cor_spf <- spear(bact_mat[idx_s,],as.matrix(mod_scores[idx_s,]))
cor_ibd <- spear(bact_mat[idx_i,],as.matrix(mod_scores[idx_i,]))
delta_r <- cor_ibd$r - cor_spf$r
n_tests <- length(cor_spf$r)

melt_cor <- function(r,p,cond) {
  df <- expand.grid(species=rownames(r), module=colnames(r),
                    stringsAsFactors=FALSE)
  df$r <- as.vector(r); df$p_adj <- as.vector(p)*n_tests; df$cond <- cond; df
}
all_cor <- rbind(melt_cor(cor_spf$r,cor_spf$p,"SPF"),
                 melt_cor(cor_ibd$r,cor_ibd$p,"IBD"))

# Top significant hits — prioritise IBD for main figure (DSS colitis focus)
top_ibd <- all_cor %>% filter(cond=="IBD", p_adj<0.05) %>%
  arrange(desc(abs(r))) %>% slice_head(n=6)
top_all <- all_cor %>% filter(p_adj<0.05) %>%
  arrange(desc(abs(r))) %>% group_by(cond) %>% slice_head(n=6) %>% ungroup()

cat("Setup complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SHARED THEME
# ─────────────────────────────────────────────────────────────────────────────
pub_theme <- theme_bw(base_size=9) +
  theme(plot.title       = element_text(size=9, face="bold", hjust=0.5),
        plot.tag         = element_text(size=11, face="bold"),
        legend.key.height= unit(0.4,"cm"),
        legend.key.width = unit(0.22,"cm"),
        legend.text      = element_text(size=7),
        legend.title     = element_text(size=7.5, face="bold"))

spatial_void <- theme_void(base_size=8) +
  theme(plot.title       = element_text(size=8, face="bold", hjust=0.5),
        plot.tag         = element_text(size=11, face="bold"),
        strip.text       = element_text(size=8, face="bold"),
        legend.key.height= unit(0.4,"cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.text      = element_text(size=6.5),
        legend.title     = element_text(size=7, face="bold"))

# ═════════════════════════════════════════════════════════════════════════════
# ── FIGURE X ─────────────────────────────────────────────────────────────────
# Panel A: Correlation heatmap (geom_tile) — SPF and IBD, side by side
# Panel B: Bivariate spatial maps for top 2 IBD hits (bact | host | overlay)
# Panel C: Quartile boxplots — module score by bacterial abundance, SPF vs IBD
# Panel D: Top host DEGs near top bacterium (bars)
# ═════════════════════════════════════════════════════════════════════════════

# ── Panel A: Correlation tile heatmap ────────────────────────────────────────
heat_df <- all_cor %>%
  mutate(stars = case_when(p_adj<0.001 ~ "***",
                            p_adj<0.01  ~ "**",
                            p_adj<0.05  ~ "*",
                            TRUE        ~ ""),
         cond = factor(cond, levels=c("SPF","IBD")))

r_lim <- max(abs(range(heat_df$r, na.rm=TRUE)), 0.15)

pA <- ggplot(heat_df, aes(x=module, y=species, fill=r)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=stars), size=2.8, color="black", vjust=0.75) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0, limits=c(-r_lim, r_lim), oob=squish,
                       name="Spearman r") +
  facet_wrap(~cond, ncol=2) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  pub_theme +
  theme(axis.text.x    = element_text(size=7, angle=40, hjust=1, color="black"),
        axis.text.y    = element_text(size=7, face="bold.italic", color="black"),
        strip.text     = element_text(size=8.5, face="bold"),
        legend.position= "right",
        panel.border   = element_rect(color="grey60", fill=NA)) +
  labs(x=NULL, y=NULL,
       title="Bacterial–host module spatial correlations")

# ── Panel B: Bivariate spatial maps (top 2 IBD hits) ─────────────────────────
make_bivar_row <- function(sp_name, mod_name, r_val, condition) {
  idx    <- all_df$condition == condition
  b_vals <- bact_mat[idx, sp_name]
  m_vals <- mod_scores[[mod_name]][idx]
  xc     <- all_df$x_coord[idx]
  yc     <- -all_df$y_coord[idx]

  s01 <- function(v) (v-min(v,na.rm=TRUE))/(diff(range(v,na.rm=TRUE))+1e-9)
  biv <- rgb(s01(b_vals), 0, s01(m_vals))

  df <- data.frame(x=xc, y=yc, b=b_vals, m=m_vals, biv=biv)

  bq <- quantile(b_vals, c(0.02,0.98), na.rm=TRUE)
  mq <- quantile(m_vals, c(0.02,0.98), na.rm=TRUE)

  p1 <- ggplot(df,aes(x,y,color=b)) +
    geom_point(size=0.35,alpha=0.85) +
    scale_color_gradientn(colours=c("grey92","#FDD49E","#FC8D59","#D7301F","#7F0000"),
                          limits=bq, oob=squish, name="log(UMI+1)") +
    coord_fixed() + spatial_void +
    labs(title=sp_name) +
    theme(plot.title=element_text(face="bold.italic"))

  p2 <- ggplot(df,aes(x,y,color=m)) +
    geom_point(size=0.35,alpha=0.85) +
    scale_color_gradientn(colours=c("grey92","#C6DBEF","#4292C6","#08306B"),
                          limits=mq, oob=squish, name="Score") +
    coord_fixed() + spatial_void +
    labs(title=paste0(mod_name," module"))

  p3 <- ggplot(df,aes(x,y)) +
    geom_point(color=df$biv, size=0.35, alpha=0.9) +
    coord_fixed() + spatial_void +
    labs(title="Co-localization") +
    annotate("text",x=-Inf,y=-Inf,hjust=-0.05,vjust=-0.4,
             label="purple = co-enriched",size=2,color="grey35")

  wrap_plots(p1,p2,p3,nrow=1) +
    plot_annotation(subtitle=paste0("[",condition,"]  r = ",round(r_val,3)),
                    theme=theme(plot.subtitle=element_text(size=7.5,hjust=0.5,
                                                           face="bold",color="grey30")))
}

# Pick top 2 IBD hits; fall back to top overall if IBD insufficient
h1 <- if(nrow(top_ibd)>=1) top_ibd[1,] else top_all[1,]
h2 <- if(nrow(top_ibd)>=2) top_ibd[2,] else top_all[2,]

pB1 <- make_bivar_row(h1$species, h1$module, h1$r, h1$cond)
pB2 <- make_bivar_row(h2$species, h2$module, h2$r, h2$cond)
pB  <- pB1 / pB2

# ── Panel C: Quartile boxplots for top 2 hits (SPF + IBD overlaid) ───────────
quart_rows <- lapply(list(h1,h2), function(h) {
  do.call(rbind, lapply(c("SPF","IBD"), function(co) {
    idx <- all_df$condition==co
    b   <- bact_mat[idx, h$species]
    m   <- mod_scores[[h$module]][idx]
    qs  <- quantile(b,0:4/4,na.rm=TRUE)
    if(length(unique(qs))<5)
      qs <- c(qs[1],quantile(b[b>0],1:3/3,na.rm=TRUE),qs[5])
    data.frame(
      Q       = cut(b,breaks=unique(qs),include.lowest=TRUE,
                    labels=paste0("Q",seq_len(length(unique(qs))-1))),
      score   = m,
      cond    = co,
      pair    = paste0(h$species,"\n× ",h$module)
    )
  }))
})
quart_df <- do.call(rbind,quart_rows)
quart_df$cond <- factor(quart_df$cond, levels=c("SPF","IBD"))
quart_df <- quart_df[!is.na(quart_df$Q),]

pC <- ggplot(quart_df, aes(Q, score, fill=cond)) +
  geom_boxplot(outlier.size=0.2, linewidth=0.35, width=0.65,
               position=position_dodge(0.75)) +
  scale_fill_manual(values=c(SPF="#4DAF4A",IBD="#E41A1C"), name=NULL) +
  facet_wrap(~pair, scales="free_y", ncol=2) +
  pub_theme +
  theme(strip.text=element_text(size=7.5,face="bold.italic"),
        legend.position="top",
        axis.text.x=element_text(size=7)) +
  labs(x="Bacterial abundance quartile",
       y="Host module score",
       title="Module score across bacterial abundance quartiles")

# ── Panel D: Top host DEGs near top bacterium ─────────────────────────────────
deg_file <- "host_bacteria_conditional_DEG.tsv"

make_pD <- function(h) {
  if (file.exists(deg_file)) {
    cdeg <- read.table(deg_file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    sub  <- cdeg[cdeg$species==h$species & cdeg$condition==h$cond, ]
  } else {
    gvar  <- apply(host_norm,2,var,na.rm=TRUE)
    tg    <- names(sort(gvar,decreasing=TRUE))[1:min(2000,ncol(host_norm))]
    idx   <- all_df$condition==h$cond
    prs   <- all_df[[h$species]][idx] > 0
    if(sum(prs)<10||sum(!prs)<10) return(NULL)
    hh    <- host_norm[idx,tg,drop=FALSE]
    res   <- do.call(rbind,lapply(tg,function(g){
      w   <- suppressWarnings(wilcox.test(hh[prs,g],hh[!prs,g],exact=FALSE))
      data.frame(gene=g,species=h$species,condition=h$cond,
                 log2FC=log2((mean(hh[prs,g],na.rm=TRUE)+1e-6)/
                             (mean(hh[!prs,g],na.rm=TRUE)+1e-6)),
                 p_val=w$p.value,stringsAsFactors=FALSE)
    }))
    res$p_adj <- p.adjust(res$p_val,"BH")
    sub <- res[order(res$p_adj),]
  }

  plot_genes <- sub %>%
    filter(p_adj<0.05, abs(log2FC)>0.3) %>%
    slice_max(order_by=abs(log2FC), n=15) %>%
    mutate(dir=ifelse(log2FC>0,"Higher","Lower"),
           gene=reorder(gene,log2FC))

  if(nrow(plot_genes)==0) return(NULL)

  ggplot(plot_genes, aes(gene, log2FC, fill=dir)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=0, linewidth=0.3) +
    scale_fill_manual(values=c(Higher="#B2182B",Lower="#2166AC"), name=NULL) +
    coord_flip() +
    pub_theme +
    theme(legend.position="top",
          axis.text.y=element_text(size=7, face="italic")) +
    labs(x=NULL,
         y="log2FC (bacteria-present vs absent)",
         title=paste0("Host DEGs near ",h$species,"\n[",h$cond,"]"))
}

pD <- make_pD(h1)
if(is.null(pD)) pD <- ggplot()+theme_void()+
  labs(title="DEG data not available — run script 4 first")

# ── Assemble Figure X ─────────────────────────────────────────────────────────
figX <- (pA | pB) / (pC | pD) +
  plot_layout(heights=c(1.6, 1), widths=c(1,1)) +
  plot_annotation(
    title  = "Host–bacteria spatial co-localization in DSS colitis",
    caption= paste0("Spearman r across spatial spots; Bonferroni correction; ",
                    "* p<0.05  ** p<0.01  *** p<0.001\n",
                    "Bivariate maps: red=bacteria, blue=host module, purple=co-enriched"),
    tag_levels = "A",
    theme = theme(
      plot.title   = element_text(size=12,face="bold",hjust=0.5),
      plot.caption = element_text(size=7,color="grey40",hjust=0)
    )
  )

ggsave("Figure_X_host_bacteria_colocalization.pdf", figX,
       width=18, height=22, units="cm", useDingbats=FALSE)
cat("Saved: Figure_X_host_bacteria_colocalization.pdf\n")

# ═════════════════════════════════════════════════════════════════════════════
# ── SUPPLEMENTARY FIGURE Y ────────────────────────────────────────────────────
# Page 1 (S-A/B): Delta correlation heatmap + ranked lollipop
# Page 2 (S-C):   All host module spatial maps (multi-panel grid)
# Page 3 (S-D/E): Spot heatmap ordered by Inflammation + cluster spatial map
# Page 4 (S-F/G): Diversity vs modules + tissue zone profiles
# ═════════════════════════════════════════════════════════════════════════════
pdf("Figure_Y_supplementary.pdf", width=18/2.54, height=24/2.54,
    useDingbats=FALSE, onefile=TRUE)

# ── Page 1: Delta heatmap (S-A) + SPF/IBD ranked lollipop (S-B) ──────────────
d_lim <- max(abs(range(delta_r,na.rm=TRUE)),0.05)

delta_df <- expand.grid(species=rownames(delta_r), module=colnames(delta_r),
                        stringsAsFactors=FALSE)
delta_df$delta <- as.vector(delta_r)
# mark pairs sig in at least one condition
delta_df$any_sig <- as.vector(
  (cor_spf$p * n_tests < 0.05) | (cor_ibd$p * n_tests < 0.05))

pSA <- ggplot(delta_df, aes(module, species, fill=delta)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(data=delta_df[delta_df$any_sig,],
            aes(label=ifelse(delta>0,"▲","▼")),
            size=2.5, color="grey20") +
  scale_fill_gradient2(low="#4DAC26",mid="white",high="#D01C8B",
                       midpoint=0,limits=c(-d_lim,d_lim),oob=squish,
                       name="Δr\n(IBD−SPF)") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  pub_theme +
  theme(axis.text.x=element_text(size=7,angle=40,hjust=1,color="black"),
        axis.text.y=element_text(size=7,face="bold.italic",color="black"),
        panel.border=element_rect(color="grey60",fill=NA),
        legend.position="right") +
  labs(x=NULL,y=NULL,
       title="S-A  Disease-driven shift in host–bacteria associations\n(▲ = strengthened in IBD; ▼ = weakened)")

# Ranked lollipop for IBD
loll_ibd <- all_cor %>% filter(cond=="IBD") %>%
  mutate(pair =paste0(species,"×",module),
         color=case_when(p_adj<0.05&r>0~"Pos*",p_adj<0.05&r<0~"Neg*",TRUE~"NS"))

pSB <- ggplot(loll_ibd, aes(reorder(pair,r),r,color=color)) +
  geom_hline(yintercept=0,linewidth=0.35,color="grey50") +
  geom_segment(aes(xend=pair,yend=0),linewidth=0.45,alpha=0.7) +
  geom_point(aes(size=color!="NS")) +
  scale_color_manual(values=c("Pos*"="#B2182B","Neg*"="#2166AC","NS"="grey65"),name=NULL) +
  scale_size_manual(values=c("TRUE"=2,"FALSE"=1),guide="none") +
  pub_theme +
  theme(axis.text.x=element_text(size=4.5,angle=60,hjust=1),
        legend.position="top") +
  labs(x=NULL,y="Spearman r",
       title="S-B  All bacterial × host module correlations (IBD, ranked)")

pg1 <- pSA / pSB + plot_layout(heights=c(1.2,1))
print(pg1)

# ── Page 2: Host module spatial maps ─────────────────────────────────────────
mod_spatial <- cbind(all_df[,c("condition","x_coord","y_coord")],mod_scores)
mod_long    <- pivot_longer(mod_spatial,cols=names(modules),
                            names_to="module",values_to="score")
mod_long$condition <- factor(mod_long$condition,levels=c("SPF","IBD"))

# Build individual module panels (2 conditions per module, returned as list)
mod_panels <- lapply(names(modules), function(md) {
  sub  <- mod_long[mod_long$module==md,]
  qlim <- quantile(sub$score,c(0.02,0.98),na.rm=TRUE)
  ggplot(sub,aes(x_coord,-y_coord,color=score)) +
    geom_point(size=0.3,alpha=0.8) +
    scale_color_gradientn(colours=c("grey92","#C6DBEF","#4292C6","#08306B"),
                          limits=qlim,oob=squish,name="Score") +
    facet_wrap(~condition,nrow=1) +
    coord_fixed() + spatial_void +
    labs(title=md) +
    theme(legend.position="right",plot.margin=margin(2,2,2,2))
})

pg2 <- wrap_plots(mod_panels, ncol=2) +
  plot_annotation(
    title="S-C  Host gene module spatial expression (SPF | IBD)",
    tag_levels=list(letters[seq_along(mod_panels)]),
    theme=theme(plot.title=element_text(size=10,face="bold",hjust=0.5)))
print(pg2)

# ── Page 3: Spot heatmap (S-D) + Cluster map (S-E) ───────────────────────────
# Spot heatmap — done via pheatmap directly (page 3 of the PDF)
if ("Inflammation" %in% names(modules)) {
  for (cond in c("SPF","IBD")) {
    idx <- which(all_df$condition==cond)
    ord <- order(mod_scores$Inflammation[idx],decreasing=TRUE)
    bh  <- t(bact_mat[idx[ord],])
    if (ncol(bh)>500) {
      keep <- round(seq(1,ncol(bh),length.out=500))
      bh   <- bh[,keep,drop=FALSE]; ord2 <- ord[keep]
    } else ord2 <- ord

    ann <- data.frame(
      Inflammation=mod_scores$Inflammation[idx[ord2]],
      row.names=seq_len(length(ord2)))
    if("Barrier" %in% names(modules))
      ann$Barrier <- mod_scores$Barrier[idx[ord2]]
    if("AntInflammatory" %in% names(modules))
      ann$AntInflammatory <- mod_scores$AntInflammatory[idx[ord2]]
    colnames(bh) <- seq_len(ncol(bh))

    ann_cols_all <- list(
      Inflammation   =colorRampPalette(c("grey95","#B2182B"))(50),
      Barrier        =colorRampPalette(c("grey95","#2166AC"))(50),
      AntInflammatory=colorRampPalette(c("grey95","#1A9641"))(50)
    )
    ann_cols <- ann_cols_all[names(ann_cols_all) %in% colnames(ann)]

    pheatmap(bh,
      cluster_rows=TRUE, cluster_cols=FALSE, show_colnames=FALSE,
      annotation_col=ann, annotation_colors=ann_cols,
      color=colorRampPalette(c("grey97","#FDD49E","#FC8D59","#D7301F","#7F0000"))(100),
      main=paste0("S-D  ",cond,
                  ": Bacterial species across spots\n",
                  "(ordered left→right by Inflammation score)"),
      fontsize=8, fontsize_row=7, border_color=NA, silent=TRUE)
  }
}

# Cluster map
set.seed(42)
feat <- scale(cbind(bact_mat,as.matrix(mod_scores)))
feat[!is.finite(feat)] <- 0
k  <- 6
km <- kmeans(feat,centers=k,nstart=25,iter.max=150)
all_df$cluster <- factor(km$cluster)

clust_pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628")

pSE1 <- ggplot(all_df,aes(x_coord,-y_coord,color=cluster)) +
  geom_point(size=0.4,alpha=0.85) +
  scale_color_manual(values=clust_pal[1:k],name="Cluster") +
  facet_wrap(~condition,nrow=1) +
  coord_fixed() + spatial_void +
  labs(title=paste0("S-E  Spatial niches (k=",k,
                    " clusters, combined host+bacterial features)")) +
  theme(legend.position="right",
        plot.title=element_text(size=9,face="bold",hjust=0.5))
print(pSE1)

# Cluster profile
cpmat <- do.call(rbind,lapply(1:k,function(ci)
  colMeans(cbind(bact_mat,as.matrix(mod_scores))[km$cluster==ci,,drop=FALSE])))
rownames(cpmat) <- paste0("C",1:k)
cpmat_z <- scale(cpmat); cpmat_z[!is.finite(cpmat_z)] <- 0

pheatmap(t(cpmat_z),
  color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  main="S-E  Cluster feature profiles (z-scored)",
  fontsize=8, fontsize_row=7, border_color="grey85",
  gaps_row=length(sp_names), silent=TRUE)

# ── Page 4: Diversity (S-F) + Tissue zones (S-G) ─────────────────────────────
props   <- sweep(bact_mat+1e-9,1,rowSums(bact_mat+1e-9),"/")
shannon <- -rowSums(props*log(props),na.rm=TRUE)

div_df   <- cbind(data.frame(shannon=shannon,condition=all_df$condition),mod_scores)
div_long <- pivot_longer(div_df,cols=names(modules),
                         names_to="module",values_to="score")
div_long$condition <- factor(div_long$condition,levels=c("SPF","IBD"))

div_ann <- div_long %>%
  group_by(condition,module) %>%
  summarise(r  =round(cor(shannon,score,method="spearman",use="complete.obs"),2),
            p  =cor.test(shannon,score,method="spearman",exact=FALSE)$p.value,
            .groups="drop") %>%
  mutate(label=paste0("r=",r,ifelse(p<0.05,"*","")),
         vjust=ifelse(condition=="SPF",1.8,3.2))

pSF <- ggplot(div_long,aes(shannon,score,color=condition)) +
  geom_point(size=0.25,alpha=0.3) +
  geom_smooth(method="lm",se=TRUE,linewidth=0.7) +
  geom_text(data=div_ann,aes(label=label,color=condition,vjust=vjust),
            x=Inf,y=Inf,hjust=1.1,size=2.4,fontface="bold",show.legend=FALSE) +
  scale_color_manual(values=c(SPF="#4DAF4A",IBD="#E41A1C"),name=NULL) +
  facet_wrap(~module,scales="free_y",ncol=4) +
  pub_theme +
  theme(strip.text=element_text(size=7.5,face="bold"),
        legend.position="top") +
  labs(x="Bacterial Shannon diversity (per spot)",y="Host module score",
       title="S-F  Bacterial community diversity vs host gene programs")

# Tissue zone profiles
zone_df <- all_df[,c("condition","x_coord","y_coord")]
zone_df$zone <- NA
for (co in c("SPF","IBD")) {
  idx <- zone_df$condition==co
  yq  <- quantile(zone_df$y_coord[idx],0:4/4,na.rm=TRUE)
  zone_df$zone[idx] <- cut(zone_df$y_coord[idx],breaks=unique(yq),
                            include.lowest=TRUE,
                            labels=paste0("Z",seq_len(length(unique(yq))-1)))
}
zone_all  <- cbind(zone_df,mod_scores)
zone_long <- pivot_longer(zone_all,cols=names(modules),
                          names_to="module",values_to="score")
zone_long$condition <- factor(zone_long$condition,levels=c("SPF","IBD"))

zone_sum <- zone_long %>%
  group_by(condition,zone,module) %>%
  summarise(mean=mean(score,na.rm=TRUE),.groups="drop")

pSG <- ggplot(zone_sum,aes(zone,mean,color=condition,group=condition)) +
  geom_line(linewidth=0.7) + geom_point(size=1.8) +
  scale_color_manual(values=c(SPF="#4DAF4A",IBD="#E41A1C"),name=NULL) +
  facet_wrap(~module,scales="free_y",ncol=4) +
  pub_theme +
  theme(strip.text=element_text(size=7.5,face="bold"),
        axis.text.x=element_text(size=7),
        legend.position="top") +
  labs(x="Tissue zone (Z1 = lumen-proximal)",y="Mean module score",
       title="S-G  Host module scores across tissue depth zones")

pg4 <- pSF / pSG + plot_layout(heights=c(1,1))
print(pg4)

dev.off()
cat("Saved: Figure_Y_supplementary.pdf\n")

cat("\n═══ Publication figures complete ═══\n")
cat("Main:          Figure_X_host_bacteria_colocalization.pdf\n")
cat("Supplementary: Figure_Y_supplementary.pdf\n")
