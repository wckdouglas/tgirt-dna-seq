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

xlim <- 720
#xlim <- 420
nucleo_p <- ggplot(data = df, aes(x = distance, weights=count/100000)) +
    geom_histogram(fill='black', color = 'black',binwidth = 6, alpha=0.8) +
    scale_x_continuous(breaks = seq(-xlim,xlim,120),limits=c(-xlim,xlim)) +
    theme(text = element_text(size=30, family='Arial', face='plain')) +
    theme(axis.text.x = element_text(size=30,face='plain', family='Arial',angle=50, hjust=0.5, vjust=0.5)) +
    theme(axis.text.y = element_text(size=30,face='plain', family='Arial')) +
    labs(x = 'Inter-nucleosome distance between\ndifferent male individual analyzed by\nTGIRT-seq and ssDNA-seq\n(bp)', 
         y = expression(Peak~count~(x10^{5}))) 
label <- expression(paste('x10'^{4}))
nucleo_p <- ggdraw(nucleo_p) 
figurename <- str_c(datapath, '/predictedNucleosomeDistance.pdf')
ggsave(nucleo_p, file=figurename, width = 11,height = 10)
message('Plotted: ', figurename)


distance_probability <-df %>% 
    mutate(filtered = case_when(
                        abs(.$distance) <= 50 ~ 'one nucleosome',
                        abs(.$distance) > 50 & abs(.$distance <= 240) ~ 'two nucleosome',
                        TRUE ~ 'far away'))  %>% 
    group_by(filtered) %>%
    summarize(total = sum(count)) %>%
    ungroup() %>%
    mutate(fraction = total/sum(total))