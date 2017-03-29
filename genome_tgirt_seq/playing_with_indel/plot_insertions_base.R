#!/usr/bin/env Rscript

library(readr)
library(stringr)
library(cowplot)
library(dplyr)
library(tidyr)
library(purrr)

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/base_insertion_table'
datapath <- '/stor/work/Lambowitz/Data/archived_work/TGIRT_ERCC_project/base_insertion_table'
table_files <- list.files(path = datapath, pattern = '.tsv', full.names = T)

paste_and_split <- function(x){
    patterns <- x$patterns
    patterns <- ifelse(!is.na(patterns),patterns,'NA')
    patterns <- str_c(patterns,collapse=',')
    patterns <- str_split(patterns,',')[[1]]
    coverage <- sum(x$coverage)
    data.frame(patterns = patterns) %>%
    mutate(coverage = coverage) %>%
    group_by(patterns, coverage) %>%
    summarize(freq = n()) %>%
    ungroup() %>%
    return()
}

read_df <- function(table_file){
    read_tsv(table_file,
             col_type = 'iccicc') %>% 
    select(-strand) %>%
    group_by(run_length, mononucleotide, indel) %>%
    nest() %>%
    mutate(pattern = map(data, paste_and_split)) %>%
    unnest(pattern) %>%
    ungroup() %>%
    filter(patterns != 'NA', patterns != '') %>%
    mutate(samplename = basename(str_replace(table_file,'.tsv',''))) %>%
    return()
}

indel_base_df <- table_files%>%
    map_df(read_df)  %>%
    group_by(indel,run_length, mononucleotide, patterns,samplename) %>%
    summarize(
        coverage = sum(coverage),
        freq = sum(freq)
    ) %>%
    ungroup() %>%
    filter(indel!='deletions') %>%
    mutate(rate = freq/coverage) %>%
    tbl_df


insert_p<-ggplot(data = indel_base_df, aes(x = run_length, y = rate,
                         color = patterns, label = patterns)) +
    geom_text() +
    facet_grid(indel~mononucleotide, scale='free_y') +
    panel_border() +
    labs(color = ' ', x = 'Homopolymer length (nt)',
         y = 'Indel rate',
         shape = '') +
    theme(text = element_text(size=30,  family = 'Arial')) +
    theme(strip.text = element_text(size = 30,  family = 'Arial'))+
    theme(axis.text = element_text(size = 30, family = 'Arial'))+
    theme(legend.key.size = unit(2,'line')) +
#    scale_color_manual(values=colors, 
#                       guide= guide_legend(ncol=2)) +
#    theme(legend.position = c(0.2,0.9)) +
#    scale_x_continuous(breaks = 4:9, labels=4:9)+
    theme(legend.position = 'none') 

figurename <- str_c(datapath, '/rna_insertion_bases.pdf')
ggsave(insert_p, file =  figurename, height = 10, width = 20)
message('plotted: ', figurename)

uniqchars <- function(x) unique(strsplit(x, "")[[1]]) 
indel_base_df %>% 
    filter(indel=='insertions') %>% 
    group_by(mononucleotide,patterns) %>% 
    summarize(freq = sum(freq)) %>%
    ungroup() %>%
    group_by(mononucleotide) %>%
    do(data_frame(
        fraction = .$freq/sum(.$freq),
        freq = .$freq,
        patterns = .$patterns
    )) %>%
    filter(mononucleotide==patterns | patterns == str_c(mononucleotide,mononucleotide))
