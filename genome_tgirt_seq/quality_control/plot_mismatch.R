#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(cowplot)


read_file <- function(tablename){
    samplename <- str_replace(basename(tablename),'.xls','')
    df <- read_tsv(tablename, skip = 1) %>%
        gather(mismatch, count, -read_pos, -sum) %>%
        separate(col = mismatch, into = c('from','to'),sep = '2') %>%
        group_by(from, to) %>%
        summarize(sum = sum(sum),
                  count = sum(count)) %>%
        mutate(mismatch_rate = count/sum) %>%
        mutate(samplename = samplename) %>%
        ungroup() %>%
        tbl_df
    return (df)
}

table_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/mismatch_profiles'
tablenames <- list.files(path = table_path, 
                         pattern = '.xls', 
                         full.names = T)
df <- tablenames %>%
    map_df(read_file)