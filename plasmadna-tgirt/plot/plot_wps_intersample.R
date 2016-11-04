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
    read_tsv()

xlim=300
nucleo_p <- ggplot(data = df, aes(x = distance, y=..count..)) +
    geom_density(fill='darkblue', color = 'darkblue') +
    scale_x_continuous(breaks = seq(-xlim,xlim,50),limits=c(-xlim,xlim)) +
    theme(text = element_text(size=35, family='Arial', face='bold')) +
    theme(axis.text.x = element_text(size=35,face='plain', family='Arial',angle=50, hjust=0.5, vjust=0.5)) +
    theme(axis.text.y = element_text(size=35,face='plain', family='Arial')) +
    labs(x = 'Distance to the Nearest Nucleosome Center (bp) \n[ssDNA-seq & TGIRT-seq]', y = 'Peak Count')
#figurename <- str_c(datapath, '/predictedNucleosomeDistance.pdf')
#ggsave(p, file=figurename)
#message('Plotted: ', figurename)
