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

xlim=400
nucleo_p <- ggplot(data = df, aes(x = distance, weights=count)) +
    geom_histogram(fill='salmon', color = 'salmon',binwidth = 3) +
    scale_x_continuous(breaks = seq(-xlim,xlim,80),limits=c(-xlim,xlim)) +
    theme(text = element_text(size=35, family='Arial', face='bold')) +
    theme(axis.text.x = element_text(size=35,face='plain', family='Arial',angle=50, hjust=0.5, vjust=0.5)) +
    theme(axis.text.y = element_text(size=35,face='plain', family='Arial')) +
    labs(x = 'Distance to the Nearest Nucleosome Center (bp) \n[ssDNA-seq & TGIRT-seq]', y = 'Peak Count')
#figurename <- str_c(datapath, '/predictedNucleosomeDistance.pdf')
#ggsave(p, file=figurename)
#message('Plotted: ', figurename)
