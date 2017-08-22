library(data.table)
library(optparse)
library(dplyr)
options(stringsAsFactors=FALSE)

option_list <- list(
  make_option(c("--counts"), type="character", default=NULL,
              help="Path to the counts_perind file from leafcutter.", metavar = "FILE"),
  make_option(c("--genotypes"), type="character", default=NULL,
              help="Path to a file with sample names and genotype counts, one per line", metavar = "FILE"),
  make_option(c("--cluster"), type="character", default=5,
              help="Name of the cluster to get counts for, e.g. clu_49269", metavar = "STRING"),
  make_option(c("--out"), type="character", default="",
              help="Base path for output files", metavar = "STRING")
)
opt <- parse_args(OptionParser(option_list=option_list))

write(paste0("Counts file: ", opt$counts), stderr())
write(paste0("Genotypes file: ", opt$genotypes), stderr())
write(paste0("cluster: ", opt$cluster), stderr())

genotypes = read.delim(opt$genotypes)

counts_perind = fread(paste0("gunzip -c ", opt$counts), header=T, showProgress=F, verbose=F, data.table=F)
names = colnames(counts_perind)[2:ncol(counts_perind)]
samplenames = sapply(strsplit(names, "[_|.]", perl=T), function(x) x[[1]])
colnames(counts_perind) = c("cluster", samplenames)

chroms = sapply(strsplit(counts_perind$cluster, ":", fixed=T), function(x) x[[1]])
starts = sapply(strsplit(counts_perind$cluster, ":", fixed=T), function(x) x[[2]])
ends = sapply(strsplit(counts_perind$cluster, ":", fixed=T), function(x) x[[3]])
clusters = sapply(strsplit(counts_perind$cluster, ":", fixed=T), function(x) x[[4]])

counts.clust = cbind(cluster=clusters, chr=chroms, start=starts, end=ends, counts_perind[, 2:ncol(counts_perind)]) %>%
  dplyr::filter(cluster==opt$cluster)

counts_perind$clustername = sapply(strsplit(counts_perind$cluster, ":", fixed=T), function(x) x[[4]])
counts_perind = counts_perind %>% dplyr::filter(clustername==opt$cluster) %>% dplyr::select(-clustername)


# Divide counts by the sum counts_perind each column
ratios.clust = sweep(counts.clust[5:ncol(counts.clust)], 2, apply(counts.clust[5:ncol(counts.clust)], 2, sum), '/')

counts.clust$gt0 = apply(counts.clust[, genotypes$sample_id[genotypes$count == 0]], MARGIN=1, sum)
counts.clust$gt1 = apply(counts.clust[, genotypes$sample_id[genotypes$count == 1]], MARGIN=1, sum)
counts.clust$gt2 = apply(counts.clust[, genotypes$sample_id[genotypes$count == 2]], MARGIN=1, sum)

counts.clust$gt0_ratio = apply(ratios.clust[, genotypes$sample_id[genotypes$count == 0]], 1, FUN=mean)
counts.clust$gt1_ratio = apply(ratios.clust[, genotypes$sample_id[genotypes$count == 1]], 1, FUN=mean)
counts.clust$gt2_ratio = apply(ratios.clust[, genotypes$sample_id[genotypes$count == 2]], 1, FUN=mean)

stderr = function(x) { sd(x) / sqrt(length(x)) }
counts.clust$gt0_ratio_se = apply(ratios.clust[, genotypes$sample_id[genotypes$count == 0]], 1, FUN=stderr)
counts.clust$gt1_ratio_se = apply(ratios.clust[, genotypes$sample_id[genotypes$count == 1]], 1, FUN=stderr)
counts.clust$gt2_ratio_se = apply(ratios.clust[, genotypes$sample_id[genotypes$count == 2]], 1, FUN=stderr)

#counts.clust$gt0_ratio = counts.clust$gt0 / sum(counts.clust$gt0)
#counts.clust$gt1_ratio = counts.clust$gt1 / sum(counts.clust$gt1)
#counts.clust$gt2_ratio = counts.clust$gt2 / sum(counts.clust$gt2)

fname = ""
if (opt$out != "") {
  fname = paste0(opt$out, ".ratios.txt")
}
rownames(ratios.clust) = counts_perind$cluster
ratios.clust = as.data.frame(t(ratios.clust))
ratios.clust$genotype = 0
ratios.clust[genotypes$sample_id[genotypes$count == 1],]$genotype = 1
ratios.clust[genotypes$sample_id[genotypes$count == 2],]$genotype = 2
write.table(ratios.clust, file=fname, sep="\t", col.names=T, row.names=F, quote=F)

if (opt$out != "") {
  fname = paste0(opt$out, ".genotype_ratios.txt")
}
write.table(counts.clust[, c("cluster", "chr", "start", "end", "gt0", "gt1", "gt2", "gt0_ratio", "gt1_ratio", "gt2_ratio", "gt0_ratio_se", "gt1_ratio_se", "gt2_ratio_se")],
            file=fname, sep="\t", col.names=T, row.names=F, quote=F)

