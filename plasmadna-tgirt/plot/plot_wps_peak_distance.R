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
    tbl_df()

dist_p <- ggplot(data = df, aes(x = distance, weight=normalized_count)) +
    geom_density(size=1.4,alpha=0.3, aes(color = prep, fill=prep,
                                         y=..count..))+
#    geom_histogram(binwidth=3, alpha=0.3, 
#                   position = 'identity',
#                   aes(y=..count../sum(..count..)*100,
#                   fill=prep, color = prep))+
    labs(x = 'Distance to the nearest\nnucleosome center (bp)', y='Percentage of peaks', color=' ', fill= ' ') +
    scale_color_manual(values =c('salmon','black')) +
    scale_fill_manual(values =c('salmon','black')) +
    theme(text = element_text(size=25, face='bold')) +
    theme(legend.key.height = unit(2,'line')) +
    theme(axis.text = element_text(size=25, face='bold')) +
    theme(legend.position = c(0.7, 0.5)) +
    xlim(0,450)
source('~/R/legend_to_color.R')
dist_p <- ggdraw(coloring_legend_text(dist_p))
figurename <- str_c(datapath, 'peak_distance.pdf',sep='/') 
ggsave(dist_p,file = figurename, width=8,height=8)
message('Plotted: ', figurename)