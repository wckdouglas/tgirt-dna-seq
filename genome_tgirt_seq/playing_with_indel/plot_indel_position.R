#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)

tablename <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/indel_pos_tables/indel_position.tsv'
df <- read_tsv(tablename) %>%
    mutate(indel_type = str_to_title(indel_type)) %>%
    group_by(indel_type, read_type, samplename) %>%
    do(data_frame(
        position = .$position,
        normalized_count = .$count / sum(.$count),
        count = .$count
    )) %>%
    filter(grepl('nextera|kh|kq',samplename)) %>%
    filter(grepl('nextera|clustered',samplename)) %>%
    mutate(prep = ifelse(grepl('nextera',samplename),'Nextera','TGIRT-seq')) %>%
    ungroup() %>%
    tbl_df

p <- ggplot(data = df, aes( x = position, y = count, 
                       group = samplename, color = prep))  +
    geom_line() +
    facet_grid(read_type~indel_type) +
    labs(x= 'Position on Read', y = 'Indel count', color = ' ')+
    theme(strip.text = element_text(face='bold',size=20)) +
    theme(text = element_text(face='bold',size=20)) +
    theme(legend.position = c(0.6,0.9))
figurename <- str_c(dirname(tablename),'/indel_position.pdf')
ggsave(p, file = figurename)