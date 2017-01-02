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
    filter(grepl('^P1|^SRR',samplename)) %>%
    filter(!grepl('cluster',samplename)) %>%
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

p <- ggplot(data = df, aes(x = distance, weight=normalized_count)) +
#    geom_density(size=1.4,alpha=0.8,aes()) +
    geom_histogram(binwidth=5, alpha=0.5, position = 'identity',aes(y=..count../sum(..count..)*100,fill=prep, color = prep))+
    labs(x = 'Distance to the Nearest Peak (bp)', y='Percentage of Peaks', color=' ', fill= ' ') +
#    scale_color_manual(values =c('salmon','black')) +
#    scale_fill_manual(values =c('salmon','black')) +
    theme(text = element_text(size=20, face='bold')) +
    theme(legend.position = c(0.7, 0.5)) +
    xlim(0,600)
figurename <- str_c(datapath, 'peak_distance.pdf',sep='/') 
ggsave(p,file = figurename)
message('Plotted: ', figurename)