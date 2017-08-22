library(data.table)
library(optparse)
options(stringsAsFactors=FALSE)

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to the counts_perind file from leafcutter.", metavar = "FILE"),
  make_option(c("--min_samples_per_intron"), type="integer", default=5,
              help="Ignore introns used (i.e. at least one supporting read) in fewer than n samples", metavar = "INT"),
  make_option(c("--min_intron_fraction"), type="numeric", default=0.01,
              help="Ignore introns with fewer than this percentage of reads across all introns in the cluster", metavar = "INT"),
  make_option(c("--min_cluster_coverage"), type="integer", default=20,
              help="Require min_samples_per_cluster samples in each cluster to have at least this many reads", metavar = "INT"),
  make_option(c("--min_samples_per_cluster"), type="integer", default=10,
              help="Require this many samples to have at least min_cluster_coverage reads in the cluster", metavar = "INT"),
  make_option(c("-o", "--outputroot"), type="character", default=NULL,
              help="Root path for output files", metavar = "PATH")
)
opt <- parse_args(OptionParser(option_list=option_list))

write(paste0("file: ", opt$file), stderr())
write(paste0("min_samples_per_intron: ", opt$min_samples_per_intron), stderr())
write(paste0("min_cluster_coverage: ", opt$min_cluster_coverage), stderr())
write(paste0("min_samples_per_cluster: ", opt$min_samples_per_cluster), stderr())
write(paste0("outputroot: ", opt$outputroot, "\n"), stderr())

counts_perind <- fread(paste0("gunzip -c ", opt$file), header=T, showProgress=F, verbose=F)

loc.df = data.frame(do.call("rbind", strsplit(counts_perind$chrom, split=":")))
#colnames(loc.df) <- c("#chr", "start", "end", "cluster")
rownames(loc.df) = counts_perind$chrom
loc.df$linenum = 1:nrow(loc.df)
counts_perind$cluster = loc.df[,4]

# N is the number of samples
N = dim(counts_perind)[2]-2
counts_perind = counts_perind[,c(1,N+2,2:(N+1)),with=F]
setkey(counts_perind, cluster)

# Output the distribution of cluster sizes
clustersize.fname = paste0(opt$outputroot, ".clustersizes.txt")
write(paste0(length(unique(counts_perind$cluster)), " clusters before filtering"), clustersize.fname)
write("Initial distribution of cluster sizes:", clustersize.fname, append=T)
write.table(table(table(counts_perind$cluster)), clustersize.fname, row.names=F, quote=F, col.names=F, append=T)

cluster.list = list()
ratios.list = list()
clusterIDs = unique(counts_perind$cluster)
for (i in 1:length(clusterIDs)) {
  clusterID = clusterIDs[i]
  if (i < 10 | i %% 100 == 0) {
    write(paste0("Processing ", i, "/", length(clusterIDs), ": ", clusterID), stderr())
  }
  cluster.data   = counts_perind[clusterID]
  cluster.counts = as.data.frame(counts_perind[clusterID, 3:(N+2), with=F])
  sample.counts = colSums(cluster.counts)
  
  samples_to_use = sample.counts > 0
  if (sum(samples_to_use) <= 1 | sum(sample.counts >= opt$min_cluster_coverage) <= opt$min_samples_per_cluster) {
    write(paste0(clusterID, ": Too few samples with cluster coverage above threshold."), stderr())
    next
  }
  intronFractions = rowSums(cluster.counts[,samples_to_use]) / sum(cluster.counts[,samples_to_use])
  introns_to_use = (rowSums(cluster.counts[,samples_to_use] > 0) >= opt$min_samples_per_intron) & (intronFractions > opt$min_intron_fraction)
  if (sum(introns_to_use) == 0) {
    write(paste0(clusterID, ": All introns of cluster filtered out."), stderr())
    next
  }
  if (sum(introns_to_use) == 1) {
    write(paste0(clusterID, ": Only 1 intron not filtered, so this cluster will be skipped."), stderr())
    next
  }
  if (any(introns_to_use == FALSE)) {
    write(paste0(clusterID, ": Skipping ", sum(!introns_to_use), " introns."), stderr())
  }
  cluster.list = c(cluster.list, list(cluster.data[introns_to_use, ]))
  
  cluster.counts = cluster.counts[introns_to_use,]
  #cluster.ratios = t(t(cluster.counts) / colSums(cluster.counts))
  cluster.ratios = sweep(cluster.counts, 2, colSums(cluster.counts), '/')
  ratios.list = c(ratios.list, list(cluster.ratios))
}

clusters.to.output = as.data.frame(rbindlist(cluster.list))
ratios.to.output = as.data.frame(rbindlist(ratios.list))

write(paste0("\n", length(unique(clusters.to.output$cluster)), " clusters after filtering"), clustersize.fname, append=T)
write("Distribution of cluster sizes after filtering (note that clusters of size 1 were removed):", clustersize.fname, append=T)
write.table(table(table(clusters.to.output$cluster)), clustersize.fname, row.names=F, quote=F, col.names=F, append=T)

clusters.to.output <- subset(clusters.to.output, select=-c(cluster))

# The ratios have the same rows in the same order as the clusters.to.output,
# but we haven't added the cluster column yet
ratios.to.output <- sapply(ratios.to.output, FUN=function(x) sprintf("%.5f", x))
ratios.to.output[as.numeric(ratios.to.output) == 0] <- "0"
ratios.to.output[as.numeric(ratios.to.output) == 1] <- "1"
ratios.to.output <- cbind(chrom=clusters.to.output[,1], as.data.frame(ratios.to.output))

# Reorder remaining clusters to original line order
clusters.to.output = clusters.to.output[order(loc.df[clusters.to.output$chrom,]$linenum),]
clusters.output.fname = paste0(opt$outputroot, ".filtered_counts.txt")
write.table(clusters.to.output, clusters.output.fname, quote=F, sep="\t", row.names=F, col.names=T)

ratios.to.output = ratios.to.output[order(loc.df[ratios.to.output$chrom,]$linenum),]
ratios.output.fname = paste0(opt$outputroot, ".ratios.txt")
write.table(ratios.to.output, ratios.output.fname, quote=F, sep="\t", row.names=F, col.names=T)

