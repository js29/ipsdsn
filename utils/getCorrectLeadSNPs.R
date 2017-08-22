library(data.table)

args <- commandArgs(trailingOnly = TRUE)
leadSnps.vcf = args[1] # VCF file from calling bcftools view passing positions of lead SNPs
fastqtl.lead.output = args[2] # Output file from Rasqual with only lead SNPs included

vcflead.df <- fread(paste0("grep -v '^##' ", leadSnps.vcf), data.table = F)
fastqtl.df <- fread(fastqtl.lead.output, data.table = F)

vcf.flt.df <- vcflead.df[vcflead.df$ID %in% fastqtl.df[,6],]

write.table(vcf.flt.df, "", quote=F, sep="\t", row.names=F, col.names=T)

