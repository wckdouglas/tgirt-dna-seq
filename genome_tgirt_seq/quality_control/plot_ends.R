#!usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)


make_prep <- function(x){
    ifelse(grepl('nextera',x),'Nextera XT','TGIRT-seq')
}

read_files <- function(filename){
    sample_name = str_replace(filename,'.csv','')
    d <- datapath %>%
        str_c(filename,sep='/') %>%
        read_csv() %>%
        mutate(filename = sample_name)
    return(d)
}

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/fragment_ends'
files <- list.files(path = datapath, pattern = '.csv')
df <- files %>%
    map(read_files) %>%
    reduce(rbind) %>%
    dplyr::filter(grepl('nextera|^K12_kh|^K12_kq', filename)) %>%
    mutate(prep = make_prep(filename)) %>% 
    mutate(read_end = ifelse(read_end == "5'", 'Read 1', 'Read 2')) %>%
    mutate(read_end = factor(read_end, levels=c("Read 1","Read 2"))) %>%
    mutate(actual_positions = ifelse(read_end == "Read 2", positions-20, positions))

p <- ggplot(data = df, aes(x = actual_positions, color = prep, linetype=filename, 
                           y = base_fraction)) +
    geom_line(size = 1.3, alpha=0.6) +
    scale_linetype_manual(guide=F,values=rep(1,length(unique(df$filename))))+
    facet_grid(base~read_end, scale ='free_x') +
    labs(x = 'Relative position to fragment ends',y='Fraction of base',color=' ') +
    panel_border() +
    scale_color_manual(values = c('light sky blue','salmon'))+
    theme(strip.text.x = element_text(size = 20, face='bold')) +
    theme(strip.text.y = element_text(size = 20, face='bold', angle = 0)) +
    theme(axis.title = element_text(size = 20, face='bold')) +
    theme(axis.text = element_text(size = 18, face='bold'))  +
    theme(legend.position = c(0.65,0.45))+
    theme(legend.text = element_text(size = 18, face='bold'))+
    theme(legend.key.size=unit(8,'mm'))
figurename <- str_c(datapath , '/end_bias_plot.pdf')
ggsave(p , file = figurename, height = 8, width = 8)
message('Plotted: ', figurename)