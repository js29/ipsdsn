library(data.table)
library(dplyr)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]

if (grepl(".gz$", infile, perl=T)) {
  dt = fread(paste("gunzip -c", infile), header = F, sep=" ", data.table=T, stringsAsFactors=F)
} else {
  dt = fread(infile, header = F, sep=" ", data.table=T, stringsAsFactors=F)  
}
colnames(dt) <- c("geneid", "ntested", "MLE1", "MLE2", "dummy", "rsid", "dist", "nominal_pval", "slope", "perm_pval", "beta_pval", "chr", "pos")

dt = na.omit(dt)

# Strangely, we have some duplicates of cluster introns / rsid combinations.
# This shouldn't happen, but my best guess is that it's due to there being
# SNPs with duplicate IDs in the VCF. To get around this we just remove the
# relatively small number of duplicates.
dt = dt[!(duplicated(paste(dt$geneid, dt$rsid))),]

loc.df = data.frame(do.call("rbind", strsplit(dt$geneid, split=":"))) 
colnames(loc.df) <- c("#chr", "start", "end", "cluster")

dt = data.table(cbind(dt, cluster=loc.df$cluster))

cluster_summary = as.data.frame(dt[, .(geneid, rsid, beta_pval, minp=min(beta_pval), numtests=length(beta_pval), bonf_pval=min(beta_pval)*length(beta_pval)), by="cluster"])
rownames(cluster_summary) = paste(cluster_summary$geneid, cluster_summary$rsid)

dt$bonf_pval = cluster_summary[paste(dt$geneid, dt$rsid), "bonf_pval"]
dt$cluster_size = cluster_summary[paste(dt$geneid, dt$rsid), "numtests"]
dt = dt[order(dt$beta_pval)]
dt = dt[!duplicated(cluster),] %>% dplyr::select(-cluster)
dt = dt[order(dt$bonf_pval)]

write.table(dt, "", quote=F, col.names=T, row.names=F, sep="\t")
