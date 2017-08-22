library(data.table)
library(optparse)
options(stringsAsFactors=FALSE)

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to the counts_perind file from leafcutter.", metavar = "FILE"),
  make_option(c("-o", "--outputroot"), type="character", default=NULL,
              help="Root path for output files", metavar = "PATH")
)
opt <- parse_args(OptionParser(option_list=option_list))

write(paste0("file: ", opt$file), stderr())
write(paste0("outputroot: ", opt$outputroot, "\n"), stderr())

counts_perind <- fread(paste0("gunzip -c ", opt$file), header=T, showProgress=F, verbose=F)

loc.df = data.frame(do.call("rbind", strsplit(counts_perind$chrom, split=":")))
colnames(loc.df) <- c("#chr", "start", "end", "cluster")
loc.df$linenum = 1:nrow(loc.df)
loc.df$start <- as.integer(loc.df$start)
loc.df$end <- as.integer(loc.df$end)

intronsizes = (loc.df$end - loc.df$start)

output.fname = paste0(opt$outputroot, ".txt")
write("summary(intronsizes):", output.fname)
write.table(data.frame(t(summary(intronsizes)))[,2:3], output.fname, row.names=F, col.names=F, quote=F, append=T)

write("quantile(intronsizes, probs=seq(0.0, 1, 0.01))", output.fname, append=T)
write.table(data.frame(quantile(intronsizes, probs=seq(0.0, 1, 0.01))), output.fname, col.names=F, quote=F, append=T)

write("\nintronsizes:", output.fname, append=T)
write.table(intronsizes, output.fname, quote=F, sep="\t", row.names=F, col.names=F, append=T)

output.fname = paste0(opt$outputroot, ".pdf")
pdf(file=output.fname, width=8, height=6)

hist(intronsizes, breaks=50, main="Intron sizes")

# Trim to get central 90% of introns
intronsizes = sort(intronsizes)
N = length(intronsizes)
bot = as.integer(0.05 * N)
top = as.integer(0.95 * N)
hist(intronsizes[bot:top], breaks=50, main="Intron sizes (central 90%)")

bot = as.integer(0.25 * N)
top = as.integer(0.75 * N)
hist(intronsizes[bot:top], breaks=50, main="Intron sizes (central 50%)")

dev.off()
