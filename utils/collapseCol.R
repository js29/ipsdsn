#!/usr/bin/env Rscript

# This function finds duplicate values based on a given set of columns,
# and then collapses another column with a separator where the first
# column are duplicates.
library(data.table)
library(dplyr)
library(optparse)

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to a file with column to collapse.", metavar = "FILE"),
  make_option(c("--header"), action="store_true", default=FALSE,
              help="Specify this option if the file has a header.", metavar = ""),
  make_option(c("--dupcols"), type="character", default=NULL,
              help="Comma-separated list of column indices, where rows may have duplicated values at these columns.", metavar = "comma-separated list of column indices"),
  make_option(c("--collapsecols"), type="character", default=NULL,
              help="Index of column(s) to collapse", metavar = "comma-separated list of column indices"),
  make_option(c("--sep"), type="character", default=",",
              help="Separator for collapsed items", metavar = "collapse separator")
)
opt <- parse_args(OptionParser(option_list=option_list))

dupcols = as.integer(strsplit(opt$dupcols, split=",")[[1]])
collapsecols = as.integer(strsplit(opt$collapsecols, split=",")[[1]])

dt = fread(opt$file, header=opt$header, data.table=T)
dt$dupcol = apply(dt[,dupcols,with=FALSE], 1, FUN=function(x) paste(x, collapse=" "))

dtnew <- dt[!duplicated(dt$dupcol),] %>% dplyr::select(-dupcol)

for (collapsecol in collapsecols) {
  collapseColName = colnames(dt)[collapsecol]
  temp = dt[, paste(eval(as.name(collapseColName)), collapse=","), by=dupcol]
  dtnew[, collapseColName] <- temp$V1  
}

write.table(dtnew, "", quote=F, col.names=opt$header, sep="\t", row.names=F)
