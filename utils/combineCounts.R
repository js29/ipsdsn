library(data.table)

args <- commandArgs(trailingOnly = TRUE)
samples.fname = args[1]
output.fname = args[2]

STARDIR = "/path/to/bams/STAR/"
sample_names = read.table(samples.fname, sep ="\t", comment.char = "", stringsAsFactors = FALSE)[,1]

loadCounts <- function(sample_dir, sample_names, counts_suffix = ".counts.txt", sub_dir = TRUE){
  #Load featureCounts output into R
  matrix = c()
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], counts_suffix, sep = ""))      
    }
    print(sample_names[i])
    
    table = fread(path, header = TRUE)
    print(head(table))
    if (i == 1){
      matrix = table[,c(1,6,7)]
    }
    else{
      matrix = cbind(matrix, table[,7])
    }
  }
  colnames(matrix) = c("gene_id", "length", sample_names)
  return(as.data.frame(matrix))
}

#Import Gencode Basic counts
basic_data = loadCounts(STARDIR, sample_names, sub_dir = TRUE, counts_suffix = ".gencode_basic.counts.txt")
write.table(basic_data, output.fname, sep = "\t", quote = FALSE, row.names = FALSE)
