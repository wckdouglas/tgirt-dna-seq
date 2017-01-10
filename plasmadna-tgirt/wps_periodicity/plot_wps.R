#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(viridis)
library(cowplot)

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA'
data_path <- str_c(project_path, '/genomeWPS/periodicity_tables')
tablename <- str_c(data_path, 'cor_periodicity_table.csv', sep='/') 
df <- tablename %>%
    read_csv() 

plot_dens <- function(label_name){
    ggplot(data=df %>% filter(label == label_name),
           aes(x=sample1,y=sample2))+
        stat_density2d(aes(fill=..level..),geom='polygon')+
        viridis::scale_fill_viridis() +
        theme(panel.background  = element_rect(fill=viridis::viridis(1))) +
        facet_wrap(~label) +
#        ylim(150,230)+
#        xlim(150,230) +
        labs(x='Nucleosome Spacing', y = 'Nucleosome Spacing', fill = 'Fraction of Window Counts') +
        theme(legend.position = 'bottom') +
        theme(strip.background = element_blank()) +
        theme(legend.text = element_text(angle=90))
}
ps <- lapply(unique(df$label), plot_dens)
p <- plot_grid(plotlist = ps)
figurename <- str_c(project_path, '/figures/periodicity_cor.pdf')
ggsave(p, file=figurename)
message('Plotted: ', figurename)