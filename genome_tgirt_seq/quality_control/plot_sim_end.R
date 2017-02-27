#!/usr/bin/env Rscript


library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)


read_files <- function(filename){
    sample_name = str_replace(filename,'.csv','')
    d <- datapath %>%
        str_c(filename,sep='/') %>%
        read_csv() %>%
        select(positions,read_end,base,base_count,base_fraction) %>%
        mutate(filename = sample_name)
    return(d)
}

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/fragment_ends'
files <- list.files(path = datapath, pattern = '.csv')
df <- files %>%
    map(read_files) %>%
    reduce(rbind) %>%
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
    mutate(prep = ifelse(sim_type=='No bias', ' ', prep)) %>%
    mutate(read_end = ifelse(read_end == "5'", 'Read 1', 'Read 2')) %>%
    mutate(read_end = factor(read_end, levels=c("Read 1","Read 2"))) %>%
    mutate(actual_positions = ifelse(read_end == "Read 2", positions-20, positions)) %>%
    mutate(read_end = as.character(read_end)) %>%
    tbl_df


p <- ggplot(data = df %>% 
                  filter(grepl('E_NEB_S6_umi2id|UMI_1_S9_umi2id|.1.MarkDuplicate',filename)))+
    geom_line(aes(x=actual_positions, y = base_fraction, color=base)) +
    geom_vline(aes(xintercept = actual_positions), linetype=2, alpha=0.3, color = 'grey') +    
    facet_grid(prep+sim_type~read_end, scale='free_x')+
    labs(x = 'Position Relative to Read ends',y='Fraction of Reads',color=' ') +
    theme(text = element_text(face='bold', size=15)) +
    panel_border()
figurename <- str_c(datapath,'/sim_frag_ends.pdf')
source('~/R/legend_to_color.R')
p <- coloring_legend_text(p)
ggsave(p, file=figurename, height = 13, width = 7)
message('Plotted: ', figurename)