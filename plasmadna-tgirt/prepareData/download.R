#!/bin/env Rscript

library(SRAdb)
library(dplyr)
library(stringr)
dbpath <- '/corral-repl/utexas/Exosome-RNA-seq/NCBI'
srafile <- str_c(dbpath,'/SRAmetadb.sqlite')
con <- dbConnect(SQLite(),srafile)
data.frame(listSRAfile('SRP061633',con)) %>%
    write.table('filelist.tsv',sep='\t',quote=F,row.names=F)
