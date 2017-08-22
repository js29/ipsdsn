library(data.table)
library(dplyr)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
sqtl_introns_fname = args[1]
exons_fname = args[2]

sqtl.introns <- fread(sqtl_introns_fname, header=T)
exons <- fread(exons_fname, header=T)
setkey(exons, EnsemblGeneID)

N = dim(sqtl.introns)[1]
for (i in 1:N) {
  intron = sqtl.introns[i,]
  geneIDs = strsplit(intron$ensemble_ids[1], ",", fixed = T)[[1]]
  symbols = strsplit(intron$symbols[1], ",", fixed = T)[[1]]
  match = NA
  for (j in 1:length(geneIDs)) {
    geneExons = exons[geneIDs[j]]
    perfectMatch = any(geneExons$ExonStart == intron$end | geneExons$ExonEnd == intron$start | geneExons$ExonStart == intron$start | geneExons$ExonEnd == intron$end)
    if (!is.na(perfectMatch)) {
      if (perfectMatch) {
        # A perfect match to one end. I assume this means we have the right gene.
        match = j
        break
      }
    }
  }
  if (!is.na(match)) {
    sqtl.introns$ensemble_ids[i] <- geneIDs[match]
    sqtl.introns$symbols[i] <- symbols[match]
  }
}

write.table(sqtl.introns, "", row.names=F, col.names=T, sep="\t", quote=F)
