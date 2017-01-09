#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)


make_prep <- function(x){
    ifelse(grepl('nextera',x),'Nextera XT','TGIRT-seq')
}


data_table <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/indel_pattern_tables/indel_pattern.tsv'
df <- read_tsv(data_table) %>%
    mutate(pattern = str_to_upper(pattern)) %>%
    mutate(indel = str_to_title(indel)) %>%
    group_by(samplename, indel) %>%
    do(data_frame(
        fraction = .$counts/sum(.$counts),
        pattern = .$pattern,
        counts = .$counts
    )) %>%
    filter(!grepl('I|D',pattern)) %>%
    ungroup() %>%
    filter(grepl('nextera|kq|kh',samplename)) %>%
    filter(grepl('nextera|cluster',samplename)) %>%
    mutate(prep = make_prep(samplename)) %>%
    mutate(first_nucleotide = str_sub(pattern,1,1)) %>%
    tbl_df


dodge <- position_dodge(width=0.9)
p <- ggplot(data = df , aes(x = pattern, y = counts, 
                            group = samplename, color = prep)) +
    geom_jitter(size = 1.5, alpha=0.6) +
    facet_grid(indel~first_nucleotide, scale='free') + 
    scale_color_manual(values = c('salmon','skyblue'))+
    labs(color= ' ', x = ' ', y = 'Count (per 10X K12 coverage)')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5 )) +
    theme(legend.position = c(0.5,0.7)) +
    theme(strip.text = element_text(size=25,face='bold'))
    
    
    