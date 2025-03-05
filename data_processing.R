tmp <- readRDS("cpg_matrix.rds")   # READ the cpg-by-cell matrix
dim(tmp)                           # 19904780 sites, 568 cells
# Filter out low-quality cells
na_by_cell <- colSums(!is.na(tmp)) # check the total cpg sites for each cell
summary(na_by_cell)                # IQR: [310983, 1272152]
ftmp <- tmp[,na_by_cell>310983]
# Filter out low-quality cpg sites
na_by_cpg <- rowSums(!is.na(ftmp))  # check total cells 
summary(na_by_cpg)                 # IQR: [8, 34], mean=25
ftmp <- ftmp[na_by_cpg>25,]
###############################################################################
# data imputation (impute by row mean)
tmp <- as.matrix(ftmp) 
tmp <- t(apply(tmp, 1, function(x) {
  x[is.na(row)] <- mean(x, na.rm = TRUE)  # Replace NA with row mean
  return(x)
}))
#save(tmp, file="~/imputed_filtered_cpg_matrix.RData")

###############################################################################
# map CpG to promoter
#load("~/imputed_filtered_cpg_matrix.RData")
library(readr)
library(GenomicRanges)
library(dplyr)
library(tidyr)

gene_data <- read.table("~/Multiomics/scnmt_parsed/features/genes/Mmusculus_genes_BioMart.87.txt", header=T)
sample_stat <- read.table("~/Multiomics/scnmt_parsed/met/results/stats/sample_stats.txt", header = T)
gene_anno <- GRanges(
  seqnames = Rle(gene_data$chr),
  ranges = IRanges(start = gene_data$start, end = gene_data$end),
  strand = Rle(gene_data$strand),
  gene_id = gene_data$symbol 
)
# Define promoter regions, 2000bp upstream from the TSS
gene_anno <- promoters(gene_anno, upstream=2000, downstream=0)

cpg_coords <- do.call(rbind, strsplit(rownames(tmp), "_"))
cpg_coords <- as.data.frame(cpg_coords, stringsAsFactors = FALSE)
colnames(cpg_coords) <- c("chr", "pos")
cpg_coords$pos <- as.numeric(cpg_coords$pos)

cpg_sites <- GRanges(
  seqnames = Rle(paste0("chr",cpg_coords$chr)),
  ranges = IRanges(start = cpg_coords$pos, end = cpg_coords$pos),
  strand = "*"
)
overlaps <- findOverlaps(cpg_sites, gene_anno)

# Map CpG sites to promoters
cpg_to_gene <- data.frame(
  cpg = rownames(tmp)[queryHits(overlaps)],
  gene_id = gene_anno$gene_id[subjectHits(overlaps)],
  stringsAsFactors = FALSE
)


tmp_df <- as.data.frame(tmp)
tmp_df$cpg <- rownames(tmp)
merged_tmp <- merge(cpg_to_gene, tmp_df, by = "cpg")

met_dt <- merged_tmp %>%
  pivot_longer(
    cols = -c(cpg, gene_id), 
    names_to = "sample", 
    values_to = "meth"
  ) %>%
  group_by(gene_id, sample) %>%
  summarize(
    sum_meth = sum(meth, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample, 
    values_from = sum_meth
  )

met_dt <- as.data.frame(met_dt) 
rownames(met_dt) <- met_dt$gene_id
met_dt <- met_dt[,-1]

#save(met_dt, file="~/promoter_matrix_sum.RData")

###############################################################################
# Filter out cells in samples with small sample sizes
sample_meta <- read.table("~/Multiomics/sample_metadata.txt", header=T)
pse.time <- read.table("~/Multiomics/real data/met_dt_pstime.txt", header=T,sep="\t")
pse.time <- merge(pse.time, sample_meta, by.x="sample", by.y = "id_rna")
pse.time <- pse.time[pse.time$id_met %in% colnames(met_dt),]

fmeta <- pse.time %>% 
  select(id_met,x,y,z,embryo,lineage10x_2)
samplecount <-fmeta %>% 
  group_by(embryo) %>% 
  summarize(n=n())
samplecount
sample_keep<-samplecount[samplecount$n>20,]
fmeta <- fmeta[fmeta$embryo %in% sample_keep$embryo,]
met_dt <- met_dt[, fmeta$id_met] #subset cells

save(met_dt, file="~/Multiomics/real data/met_dt_sum.RData")