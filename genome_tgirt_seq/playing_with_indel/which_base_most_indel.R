#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(purrr)
library(stringr)
library(tidyr)
library(extrafont)

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/base_indel_table'
indel_tables <- list.files(path = datapath, 
                           pattern = '.tsv',
                           full.names = T)
indel_tables <- indel_tables[grepl('clustered_fam',indel_tables)]
indel_tables <- indel_tables[grepl('30X',indel_tables)]

complementary_base <- function(b){
    return (str_replace(b, 'ACTG','TGAC'))
}


read_file <- function(indel_table){
    table_name <- basename(indel_table)
    samplename <- str_replace(table_name,'.tsv','')
    df <- read_tsv(indel_table) %>%
        gather(indel, counts, -mononucleotide, -run_length) %>%
        mutate(strand = ifelse(grepl('rev|neg',indel),'Reverse strand', 'Foward strand')) %>%
        mutate(indel = str_replace_all(indel,'negatve_|positive_|fwd_|rev_','')) %>%
        mutate(mononucleotide = ifelse(strand == 'Reverse strand',chartr('ACTG','TGAC',mononucleotide),mononucleotide)) %>%
        group_by(mononucleotide, run_length, indel) %>%
        summarize(counts = sum(counts)) %>%
        spread(indel, counts) %>%
        gather(indel, counts, -mononucleotide:-aln_count) %>%
        mutate(indel_rate = counts/aln_count) %>%
        mutate(method = ifelse(grepl('nextera',samplename),'Nextera-XT','TGIRT-seq')) %>%
        mutate(samplename = samplename) %>%
        filter(counts != 0) %>%
        tbl_df 
    return(df)
}

df <- indel_tables %>%
    map_df(read_file) %>%
    filter(method != 'Nextera-XT') %>%
    mutate(indel = str_to_title(indel)) %>%
    tbl_df

colors <- c('salmon','deepskyblue3','goldenrod4','darkgreen')
base_indel_p<-ggplot(data = df, 
       aes(color = mononucleotide, 
           x = run_length, 
           #shape = method,
           y = indel_rate)) + 
    geom_jitter(alpha=0.8, size=1) + 
    facet_grid(indel~mononucleotide, scale='free_y')+
    panel_border()+
    #viridis::scale_color_viridis() +
    labs(color = ' ', x = 'Mononucleotide runs (nt)',
         y = 'Indel rate',
         shape = '') +
    theme(text = element_text(size=30,  family = 'Arial')) +
    theme(strip.text = element_text(size = 30,  family = 'Arial'))+
    theme(axis.text = element_text(size = 30, family = 'Arial'))+
    theme(legend.key.size = unit(2,'line')) +
    scale_color_manual(values=colors, 
                       guide= guide_legend(ncol=2)) +
    theme(legend.position = c(0.2,0.9)) +
    theme(legend.position = 'none') +
    scale_x_continuous(breaks = 4:9, labels=4:9)
source('~/R/legend_to_color.R')
base_indel_p<-ggdraw(coloring_legend_text_match(base_indel_p,colors))
    
figure_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures'
figurename <- str_c(figure_path, '/base_mononucleotide.pdf')
ggsave(base_indel_p, file =  figurename, height = 7, width = 10)
message('plotted: ', figurename)
