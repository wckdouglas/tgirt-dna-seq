#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(cowplot)

read_file <- function(tablename){
    samplename <- str_replace(basename(tablename),'.tsv','')
    read_tsv(tablename) %>%
    select(-samplename) %>%
    mutate(filename = samplename) %>%
    return()
}

data_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/isize_tables'
tablenames <- list.files(path = data_path,
           pattern = '.tsv',
           full.names = T)

df <- tablenames %>%
    map_df(read_file)  %>%
    filter(grepl('UMI|NEB|sim',filename)) %>%
    filter(grepl('MarkDuplicate',filename)) %>%
    filter(grepl('sim|umi2id',filename)) %>%
    mutate(prep = case_when(
        grepl('13N|UMI',.$filename) ~ 'Covaris',
        grepl('fragmentase|NEB',.$filename) ~ 'Fragmentase'
    )) %>%
    mutate(sim_type = case_when(
        grepl('ligation', .$filename) ~ "5' end only",
        grepl('template', .$filename) ~ "3' end only",
        grepl('no_bias',.$filename) ~'No bias',
        grepl('sim',.$filename) ~ 'both ends',
        grepl('^K12',.$filename) ~'Experimental'
    )) %>%
    filter(!is.na(prep)) %>%
    mutate(prep = ifelse(sim_type=='No bias', 'No bias', prep)) %>%
    tbl_df

source('~/R/legend_to_color.R')
p <- ggplot(data=df, aes(x=isize,y=percentage*100, color = sim_type, group_filename)) +
    geom_line() +
    facet_grid(prep~.) +
    theme(text = element_text(face='bold', size=25)) +
    theme(axis.text = element_text(face='bold', size=25)) +
    panel_border() +
    theme(legend.key.height = unit(2,'line')) +
    labs(color = ' ', x= 'Fragment size', y='Percentage of fragments') 
p <- ggdraw(coloring_legend_text(p))
figure_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures'
figurename <- str_c(figure_path, '/sim_frag_size.pdf')
ggsave(p, file=figurename, height = 7, width = 10)
message('Plotted: ', figurename)

    