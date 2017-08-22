library("dplyr")
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
expr_dataset_rds = args[1]
expressed_only = args[2]

expr_dataset = readRDS(expr_dataset_rds)

chromosome.df = data.frame(chr=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"),
                           chrorder=seq(1,25))
rownames(chromosome.df) = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")

genepos = dplyr::select(expr_dataset$gene_metadata, chromosome_name, start_position, end_position, gene_id) %>% 
  dplyr::rename_("chr" = "chromosome_name", "left" = "start_position", "right" = "end_position") %>%
  dplyr::arrange(chr, left) %>%
  dplyr::distinct()

rownames(expr_dataset$design) <- expr_dataset$design$sample_id
hipsci_ids <- expr_dataset$design[expr_dataset$design$sample_id, "genotype_id"]

colnames(expr_dataset$exprs_cqn) <- hipsci_ids

exprs_sel = expr_dataset$exprs_cqn
if (!is.na(expressed_only) && expressed_only != FALSE) {
  expressed_genes = names(which(rowMeans(expr_dataset$exprs_cqn) > 0)) #Set conservative threshold to expression level
  exprs_sel = expr_dataset$exprs_cqn[expressed_genes,]
}

res = dplyr::mutate(as.data.frame(exprs_sel), gene_id = rownames(exprs_sel)) %>%
  dplyr::left_join(., genepos, by = "gene_id") %>%
  dplyr::select(chr, left, right, gene_id, everything()) %>%
  dplyr::left_join(., chromosome.df, by = "chr") %>%
  dplyr::arrange(chrorder, left) %>%
  dplyr::select(-chrorder) %>%
  dplyr::rename_("#chr" = "chr")

write.table(res, file="", quote=F, sep="\t", row.names=F)
