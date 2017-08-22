args <- commandArgs(trailingOnly=TRUE)
inputRds <- args[1]
outputRoot <- args[2]

vcf_file = readRDS(inputRds) #genotypes
vcf_file$genotypes$ID <- rownames(vcf_file$genotypes)
vcf_file$genotypes$chr <- vcf_file$snpspos$chr[vcf_file$genotypes$ID]

chrnames <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
for (chr in chrnames) {
  chr_genotypes <- vcf_file$genotypes[vcf_file$genotypes$chr == chr, !(colnames(vcf_file$genotypes) == "chr")]
  fname <- paste0(outputRoot, ".genotypes.", chr, ".txt")
  write.table(chr_genotypes, file=fname, row.names=F, col.names=T, quote=F, sep="\t")
  #fname <- paste0(outputRoot, ".snppos.", chr, ".txt")
  #write.table(vcf_file$snpspos[vcf_file$snpspos$chr == chr], file=fname, row.names=F, col.names=T, quote=F, sep="\t")
}

