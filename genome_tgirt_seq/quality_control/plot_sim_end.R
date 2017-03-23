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
    purrr::reduce(rbind) %>%
    filter(grepl('13N_K12_sim|13N_K12_sim_ligation_only|no_bias|UMI',filename)) %>%
    #filter(grepl('13N_kmer_K12_sim|13N_K12_sim_ligation_only|no_bias|UMI',filename)) %>%
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
                        grepl('sim',.$filename) ~ 'Both ends',
                        grepl('^K12',.$filename) ~'Experimental'
                        )) %>%
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
    labs(x = 'Position Relative to Read ends',y='Fraction of Reads',color=' ') +
    theme(text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(strip.text.y = element_text(face='bold', size=18)) +
    panel_border() +
    theme(legend.position = c(0.2,0.95))  +
    scale_color_discrete(guide = guide_legend(ncol=4))
figurename <- str_c(datapath,'/sim_frag_ends.pdf')
source('~/R/legend_to_color.R')
sim_end_p <- ggdraw(coloring_legend_text(sim_end_p))
ggsave(sim_end_p, file=figurename, height = 13, width = 7)
message('Plotted: ', figurename)

ddf <- df %>% 
    filter(grepl('E_NEB_S6_umi2id|UMI_1_S9_umi2id|.2.MarkDuplicate',filename)) %>%
    filter(grepl('sim.2|umi2id', filename)) %>% 
    mutate(sim_type = ifelse(grepl('Both',sim_type),'Simulation','Experimental')) %>%
    select(sim_type, base_fraction, positions, read_end,base, prep) %>%
    group_by(sim_type, positions, read_end, prep) %>%
#    summarize(bit = -sum(base_fraction * log2(base_fraction))) %>%
#    spread(sim_type, bit) %>%
    spread(sim_type, base_fraction) %>%
    mutate(residual =  -(Experimental - Simulation) )

dp <- ggplot(data = ddf, aes(x = positions, y = residual, color = base)) +
    geom_line() +
    geom_hline(yintercept = 0, color = 'grey', alpha=0.7, linetype='dashed')+
    facet_grid(prep~read_end, scale='free_x') +
    labs(y = 'Residual\n(Simulation - Experimental)', x='Position relative to read ends (nt)', color = ' ') +
    panel_border() +
    theme(axis.text = element_text(size = 25, face= 'bold')) +
    theme(axis.title = element_text(size = 25, face= 'bold')) +
    theme(strip.text = element_text(size = 25, face= 'bold'))
    

