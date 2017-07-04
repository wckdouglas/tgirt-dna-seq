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
simulation_annotation <- c('Experimental',"5' and 3' bias","5' bias","3' bias","No bias")
files <- list.files(path = datapath, pattern = '.csv')
df <- files %>%
    .[grepl('amp|UMI', .)] %>%
    .[grepl('13N_clustered_K12_sim|no_bias|UMI', .)] %>%
    map(read_files) %>%
    purrr::reduce(rbind) %>%
    #filter(grepl('13N_kmer_K12_sim|13N_K12_sim_ligation_only|no_bias|UMI',filename)) %>%
    filter(grepl('MarkDuplicate',filename)) %>%
    filter(grepl('sim|umi2id',filename)) %>%
    mutate(prep = case_when(
                        grepl('13N|UMI',.$filename) ~ 'Covaris',
                        grepl('fragmentase|NEB',.$filename) ~ 'Fragmentase'
                        )) %>%
    mutate(sim_type = case_when(
                        grepl('ligation', .$filename) ~ simulation_annotation[3],
                        grepl('template', .$filename) ~ simulation_annotation[4],
                        grepl('no_bias',.$filename) ~ simulation_annotation[5],
                        grepl('sim',.$filename) ~ simulation_annotation[2],
                        grepl('^K12',.$filename) ~ simulation_annotation[1]
                        )) %>%
    mutate(sim_type = factor(sim_type, levels = simulation_annotation )) %>%
    filter(!is.na(prep)) %>%
    mutate(prep = ifelse(sim_type=='No bias', ' ', prep)) %>%
    mutate(read_end = ifelse(read_end == "5'", 'Read 1', 'Read 2')) %>%
    mutate(read_end = factor(read_end, levels=c("Read 1","Read 2"))) %>%
    mutate(actual_positions = ifelse(read_end == "Read 2", positions-20, positions)) %>%
    mutate(read_end = as.character(read_end)) %>%
    tbl_df



colors <- c('red','green','blue','purple')
sim_end_p <- ggplot(data = df %>% 
                filter(!grepl('Fragmentase',prep)) %>%
                filter(grepl('E_NEB_S6_umi2id|UMI_1_S9_umi2id|.1.MarkDuplicate',filename)))+
    geom_line(aes(x=actual_positions, y = base_fraction, color=base)) +
    geom_vline(aes(xintercept = actual_positions), linetype=2, alpha=0.3, color = 'grey') +    
    facet_grid(sim_type~read_end, scale='free_x')+
    labs(x = 'Position relative to read ends',y='Fraction of Reads',color=' ') +
    theme(text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(panel.spacing = unit(2, "lines")) +
    theme(strip.text.y = element_text(face='plain',family='Arial', size=30)) +
    panel_border() +
    theme(legend.position = c(0.2,0.95))  +
    scale_color_discrete(guide = guide_legend(ncol=4))
figurename <- str_c(datapath,'/sim_frag_ends.pdf')
source('~/R/legend_to_color.R')
sim_end_p <- ggdraw(coloring_legend_text(sim_end_p))
ggsave(sim_end_p, file=figurename, height = 13, width = 7)
message('Plotted: ', figurename)
