#!/usr/bin/env Rscript

library(readr)
library(stringr)
library(cowplot)
library(dplyr)
library(tidyr)
library(purrr)

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/base_insertion_table'
table_files <- list.files(path = datapath, pattern = '.tsv', full.names = T)

paste_and_split <- function(d){
    x <- d[[1]]
    x <- str_c(x,collapse=',')
    res <- str_split(x,',')[[1]]
    return(res)
}

read_df <- function(table_file){
    read_tsv(table_file) %>% 
    filter(patterns != ',') %>% 
    mutate(mononucleotide = ifelse(strand == '-',chartr('ACTG','TGAC',mononucleotide),mononucleotide))  %>%
    select(-strand) %>%
    group_by(run_length, mononucleotide, indel) %>%
    nest() %>%
    mutate(pattern = map(data, paste_and_split)) %>%
    unnest(pattern) %>%
    filter(pattern!='') %>%
    mutate(samplename = basename(str_replace(table_file,'.tsv',''))) %>%
    return()
}

df <- table_files%>%
    map_df(read_df)  %>%
    group_by(pattern,indel,mononucleotide, run_length) %>%
    summarize(freq = n()) 


p<-ggplot(data = df, aes(x = run_length, y = freq, color = mononucleotide, label = pattern)) +
    geom_text() +
    facet_grid(indel~mononucleotide) +
    panel_border() +
    scale_y_log10()