library(Seurat)
library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
#library(wesanderson)
library(dplyr)
library(Matrix)
library(viridis)


dir <- "D:/microbiome_revision/compare_SPF_IBD" 
setwd(dir) 
#imported_raster=OpenImageR::readImage("MCSI0619_small.jpg")

SPF = read.table(file = 'GSM8921752_SPF_stdata.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE) 

IBD = read.table(file = 'GSM8921753_SPFIBD_stdata.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE) 

# prefixes already defined
prefixes <- c(
  "GOZ73", "CL2", "F2Y25", "K6V30", "I6J55", "Bovatus", "HMPREF2141",
  "RY69", "DW767", "LA357", "CSV91", "DWX94", "AQ990", "LUA64",
  "ENT", "ELI", "Fi14EGH31", "C4Q21", "DXB85", "FYC25", "DWX93", "CDL20"
)

# find all gene columns in mydata with these prefixes
target_genes <- grep(
  pattern = paste0("^(", paste(prefixes, collapse = "|"), ")_RS"),
  x = colnames(my_data),
  value = TRUE
)


# calculate total UMI across all spots for each gene
gene_umi_totals <- colSums(my_data[, target_genes, drop = FALSE], na.rm = TRUE)

# convert to data frame
gene_umi_totals_df <- data.frame(
  gene = names(gene_umi_totals),
  prefix = sub("_RS.*", "", names(gene_umi_totals)),
  total_UMI = as.numeric(gene_umi_totals),
  row.names = NULL
)

# sort from highest to lowest total UMI
gene_umi_totals_df <- gene_umi_totals_df[order(-gene_umi_totals_df$total_UMI), ]

# show result
print(gene_umi_totals_df)

prefix_totals <- aggregate(total_UMI ~ prefix, data = gene_umi_totals_df, sum)
prefix_totals <- prefix_totals[order(-prefix_totals$total_UMI), ]
print(prefix_totals)

write.table(
  prefix_totals,
  file = "prefix_totals.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



species_map <- c(
  "GOZ73"      = "A. muciniphila",
  "CL2"        = "A. hadrus",
  "F2Y25"      = "B. caccae",
  "K6V30"      = "B. faecis",
  "I6J55"      = "B. fragilis",
  "Bovatus"    = "B. ovatus",
  "HMPREF2141" = "B. uniformis",
  "RY69"       = "B. breve",
  "DW767"      = "B. obeum",
  "LA357"      = "C. sporogenes",
  "CSV91"      = "C. aerofaciens",
  "DWX94"      = "C. eutactus",
  "AQ990"      = "D. longicatena",
  "LUA64"      = "E. lenta",
  "ENT"        = "E. faecalis",
  "ELI"        = "E. callanderi",
  "Fi14EGH31"  = "F. intestinalis",
  "C4Q21"      = "F. prausnitzii",
  "DXB85"      = "P. merdae",
  "FYC25"      = "P. clara",
  "DWX93"      = "R. hominis",
  "CDL20"      = "R. gnavus"
)


