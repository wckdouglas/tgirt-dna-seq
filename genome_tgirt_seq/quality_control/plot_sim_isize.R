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

simulation_annotation <- c('Experimental',"5' and 3' bias","5' bias","3' bias","No bias")
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
        grepl('ligation', .$filename) ~ simulation_annotation[3],
        grepl('template', .$filename) ~ simulation_annotation[4],
        grepl('no_bias',.$filename) ~ simulation_annotation[5],
        grepl('sim',.$filename) ~ simulation_annotation[2],
        grepl('^K12',.$filename) ~ simulation_annotation[1]
    )) %>%
    filter(!is.na(prep)) %>%
    mutate(prep = ifelse(sim_type=='No bias', 'No bias', prep)) %>%
    filter(prep != 'Fragmentase') %>%
    tbl_df

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

colors <- c('black','red','goldenrod4','springgreen4','navyblue','grey72')
source('~/R/legend_to_color.R')
sim_size_p <- ggplot(data=df, aes(x=isize,y=percentage*100, 
                                  color = factor(sim_type, levels=simulation_annotation), 
                                  group=filename)) +
    geom_line() +
    #facet_grid(prep~.) +
    theme(text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(legend.key.height = unit(2,'line')) +
    theme(legend.position = c(0.6,0.6))+
    scale_color_manual(values = colors) + 
    labs(color = ' ', x= 'Fragment size', y='Percentage of fragments') 
sim_size_p <- ggdraw(coloring_legend_text_match(sim_size_p, colors))
figure_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures'
figurename <- str_c(figure_path, '/sim_frag_size.pdf')
ggsave(sim_size_p, file=figurename, height = 7, width = 10)
message('Plotted: ', figurename)

    