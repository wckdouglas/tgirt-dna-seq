#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(extrafont)
library(FBN)
loadfonts()

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures'
df <- datapath %>%
    str_c('predictedNucleosomeDistance.tsv',sep='/') %>%
    read_tsv() %>%
    group_by(distance) %>%
    summarize(count = n()) %>%
    ungroup()

xlim=720
nucleo_p <- ggplot(data = df, aes(x = distance, weights=count/10000)) +
    geom_histogram(fill='salmon', color = 'salmon',binwidth = 6) +
    scale_x_continuous(breaks = seq(-xlim,xlim,80),limits=c(-xlim,xlim)) +
    theme(text = element_text(size=30, family='Arial', face='bold')) +
    theme(axis.text.x = element_text(size=30,face='plain', family='Arial',angle=50, hjust=0.5, vjust=0.5)) +
    theme(axis.text.y = element_text(size=30,face='plain', family='Arial')) +
    labs(x = 'Distance to the Nearest Nucleosome Center (bp)\n[ssDNA-seq (ref.2) & TGIRT-seq]', y = 'Peak Count') 
label <- expression(paste('x10'^{'5'}))
nucleo_p <- ggdraw(nucleo_p) +
    draw_label(label, x = 0.11, y = 0.95, size = 25, fontface ='bold')
figurename <- str_c(datapath, '/predictedNucleosomeDistance.pdf')
ggsave(nucleo_p, file=figurename, width = 11,height = 10)
message('Plotted: ', figurename)
