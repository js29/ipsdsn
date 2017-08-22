library(data.table)
library(dplyr)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
# File with counts per individual per cluster intron, including cluster sums
leafcutter.ratios.fname = args[1]
sample.meta.fname = args[2]
outputRoot = args[3]

if (grepl(".gz$", leafcutter.ratios.fname, perl=T)) {
  counts.df = fread(paste("gunzip -c", leafcutter.ratios.fname), data.table=F, stringsAsFactors=F)
} else {
  counts.df = fread(leafcutter.ratios.fname, data.table=F, stringsAsFactors=F)  
}
#counts.df = fread("snqtl.v5.clusters_perind.ratios.gz.head", data.table=F, stringsAsFactors=F, sep=" ", header=T)

# Replace sample BAM file names in header with HIPSCI genotype IDs
sample.meta = read.delim(sample.meta.fname, header=T)
leafNames = colnames(counts.df)[-1]
sampleIDs = sapply(leafNames, FUN=function(y) strsplit(y, "\\.")[[1]][1])
rownames(sample.meta) = sample.meta$sample_id
colnames(counts.df)[-1] = sample.meta[sampleIDs,]$genotype_id

counts.mat.scaled = t(as.matrix(counts.df[-1], mode="integer")) %>% scale() %>% t()

quantiles.normal <- function(v) {
  ranking <- rank(v, ties.method="average")
  n <- length(ranking)
  u <- qnorm(seq(from=1/(n+1), to=n/(n+1), by=1/(n+1)))
  w <- u[ranking]
  return(w)
}
quantiles.normal.ignoreNA = function(v) {
  v[!is.na(v)] = quantiles.normal(v[!is.na(v)])
  v
}

counts.mat.scaled.qnorm = apply(counts.mat.scaled, MARGIN=2, quantiles.normal.ignoreNA)
  
loc.df = data.frame(do.call("rbind", strsplit(counts.df$chrom, split=":"))) 
colnames(loc.df) <- c("#chr", "start", "end", "cluster")
loc.df$start <- as.integer(loc.df$start)
loc.df$end <- as.integer(loc.df$end)

# FastQTL uses the START coordinate to define the window for testing SNPs.
# But we want it to use the middle of the intron, so we need to change the
# start coord to be the middle.
loc.df$start <- loc.df$start + floor((loc.df$end - loc.df$start) / 2)
loc.df$end <- loc.df$start + 1

loc.df$cluster = counts.df$chrom
counts.normalized.df <- cbind(loc.df, counts.mat.scaled.qnorm)

rownames(counts.mat.scaled.qnorm) <- loc.df$cluster

is.num <- sapply(counts.normalized.df, is.numeric)
counts.normalized.df[is.num] <- lapply(counts.normalized.df[is.num], round, 5)
write.table(counts.normalized.df, paste0(outputRoot, ".normalized.txt"), sep="\t", row.names=F, col.names=T, quote=F)

# Also write out expression PCs
counts.mat.scaled.qnorm.nona = na.omit(counts.mat.scaled.qnorm)
pca = prcomp(t(counts.mat.scaled.qnorm.nona))
pca_df = as.data.frame(t(pca$x)) %>%
  dplyr::mutate(PC = colnames(pca$x)) %>%
  dplyr::select(PC, everything())

write.table(pca_df, file=paste0(outputRoot, ".pcs.txt"), row.names=F, col.names=T, quote=F, sep="\t")

