library("dplyr")
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
expr_dataset_rds = args[1]
expressed_only = args[2]
useHipsciGT = TRUE

expr_dataset = readRDS(expr_dataset_rds)

if (useHipsciGT) {
  rownames(expr_dataset$design) <- expr_dataset$design$sample_id
  hipsci_ids <- expr_dataset$design[expr_dataset$design$sample_id, "genotype_id"]
  colnames(expr_dataset$exprs_cqn) <- hipsci_ids
}

exprs_sel = expr_dataset$exprs_cqn
if (!is.na(expressed_only) && expressed_only != FALSE) {
  expressed_genes = names(which(rowMeans(expr_dataset$exprs_cqn) > 0)) #Set conservative threshold to expression level
  exprs_sel = expr_dataset$exprs_cqn[expressed_genes,]
}

pca = prcomp(t(exprs_sel))
pca_df = as.data.frame(t(pca$x)) %>%
  dplyr::mutate(PC = colnames(pca$x)) %>%
  dplyr::select(PC, everything())

write.table(pca_df, file="", row.names=F, col.names=T, quote=F, sep="\t")
