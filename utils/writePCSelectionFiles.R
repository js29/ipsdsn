library("dplyr")
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
numPCs = as.integer(args[1])
outputRoot = args[2]

pcsVec = paste0(rep("PC", numPCs), seq(1, numPCs, 1))

for (i in 1:numPCs) {
  fname = paste0(outputRoot, ".", i, ".txt")
  write.table(pcsVec[1:i], file=fname, row.names=F, col.names=F, quote=F)
}
