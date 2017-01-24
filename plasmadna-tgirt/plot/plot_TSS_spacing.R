#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)
library(cowplot)

make_cor_df <- function(filename, datapath){
    periodogram <- datapath %>%
        str_c(filename , sep='/') %>%
        read_tsv() %>%
        mutate(samplename = str_replace(filename,'.bed','')) %>%
        tbl_df
    return(periodogram)
}

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/tss_periodicity'
files <- list.files(path = datapath, pattern='.bed', full.names = F)
df <- files %>%
    mclapply(.,make_cor_df, datapath, mc.cores=12) %>%
    reduce(rbind) %>%
    spread(samplename, intensity) %>%
    tbl_df
