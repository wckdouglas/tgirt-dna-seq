#!/usr/bin/env Rscript

library(readr)
library(purrr)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(cowplot)

indel_index_table <- '/stor/work/Lambowitz/ref/Ecoli/k12_mg1655_repeat_index.tsv'
indel_index_df <- read_tsv(indel_index_table)

read_indel_table <- function(tablename){
    samplename <- basename(str_replace(tablename,'.tsv',''))
    tablename %>%
    read_tsv() %>%
    mutate(samplename = samplename) %>%
    tbl_df %>%
    return()
}


indel_table_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/indel_table'
df <- list.files(path = indel_table_path, pattern = '.tsv', full.names = T) %>%
    map_df(read_indel_table) %>%
    inner_join(indel_index_df)%>% 
    select(grep('samplename|num_D|num_I|index',names(.))) %>%
    filter(!grepl('MiSeq|Ecoli',samplename)) %>%
    mutate(prep = ifelse(grepl('nextera',samplename), 'Nextera XT','TGIRT-seq')) %>%
    mutate(indel_index = negative_index + positive_index ) %>%
    mutate(number_of_indel = num_D + num_I) %>%
    group_by(prep, indel_index) %>%
    summarize(number_of_indel = sum(number_of_indel),
              count = n()
    ) %>%
    mutate(normalized_indel = number_of_indel / count) %>%
    tbl_df
    
form <- y ~ poly(x, 3, raw = TRUE)
p<-ggplot(data = df, aes(x = indel_index, y = normalized_indel, color = prep))+
    geom_smooth(method='lm', formula =form) +
    geom_point() +
    stat_poly_eq(aes(label = ..eq.label..), formula = form, 
                 parse = TRUE)