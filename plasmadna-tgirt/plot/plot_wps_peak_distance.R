#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(tibble)


datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures'
df <- datapath %>%
    str_c('peakDistance.tsv',sep='/') %>%
    read_tsv()  %>%
    filter(grepl('umi|0052$',samplename)) %>%
    filter(grepl('P1022|0052$',samplename)) %>%
    mutate(prep = ifelse(grepl('SRR',samplename),'ssDNA-seq','TGIRT-seq')) %>%
    filter(distance < 450) %>%
    group_by(prep,distance) %>% 
    summarize(nucleosome_count = sum(nucleosome_count)) %>%
    group_by(prep) %>%
    do(data_frame(
        normalized_count = .$nucleosome_count / sum(.$nucleosome_count),
        distance = .$distance,
        nucleosome_count = .$nucleosome_count
    ))  %>%
    ungroup() %>%
    mutate(prep = factor(prep, level = rev(unique(prep)))) %>%
    tbl_df()

dist_p <- ggplot(data = df, aes(x = distance, weight=normalized_count)) +
    geom_density(size=1.4,alpha=0.3, linetype=1, aes(color = prep,# fill=prep,
                                         y=..count..), show.legend =F)+
    stat_density(aes(x=distance, color=prep, weight=normalized_count),
                 geom="line",position="identity", size = 0) +
#    geom_histogram(binwidth=3, alpha=0.3, 
#                   position = 'identity',
#                   aes(y=..count../sum(..count..)*100,
#                   fill=prep, color = prep))+
    labs(x = 'Inter-nucleosome distance in cfDNA of two individuals (bp)', y='% Peaks', color=' ', fill= ' ') +
    scale_color_manual(values =c('black','salmon')) +
    theme(text = element_text(size=30, face='plain', family = 'Arial')) +
    theme(legend.key.height = unit(2,'line')) +
    theme(axis.text = element_text(size=30, face='plain', family = 'Arial')) +
    theme(legend.position = c(0.7, 0.5)) +
    xlim(0,450)+
    guides(color = guide_legend(override.aes=list(size=1)))
source('~/R/legend_to_color.R')
dist_p <- ggdraw(coloring_legend_text(dist_p))
figurename <- str_c(datapath, 'peak_distance.pdf',sep='/') 
ggsave(dist_p,file = figurename, width=8,height=8)
message('Plotted: ', figurename)