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

xlim=720
nucleo_p <- ggplot(data = df, aes(x = distance, weights=count/10000)) +
    geom_histogram(fill='salmon', color = 'salmon',binwidth = 6) +
    scale_x_continuous(breaks = seq(-xlim,xlim,80),limits=c(-xlim,xlim)) +
    theme(text = element_text(size=30, family='Arial', face='bold')) +
    theme(axis.text.x = element_text(size=30,face='plain', family='Arial',angle=50, hjust=0.5, vjust=0.5)) +
    theme(axis.text.y = element_text(size=30,face='plain', family='Arial')) +
    labs(x = 'Difference in distance\nbetween nucleosome centers(bp)\n[ssDNA-seq (ref.2) vs TGIRT-seq]', 
         y = 'Peak count') 
label <- expression(paste('x10'^{4}))
nucleo_p <- ggdraw(nucleo_p) +
    draw_label(label, x = 0.17, y = 0.98, size = 25, fontface ='bold')
figurename <- str_c(datapath, '/predictedNucleosomeDistance.pdf')
ggsave(nucleo_p, file=figurename, width = 11,height = 10)
message('Plotted: ', figurename)


df %>% 
    mutate(filtered = ifelse(distance < 240 & distance > -240, 'in','out'))  %>% 
    mutate(filtered = ifelse(distance ==0,'0', filtered))%>%
    group_by(filtered) %>%
    summarize(total = sum(count)) %>%
    ungroup() %>%
    mutate(fraction = total/sum(total))