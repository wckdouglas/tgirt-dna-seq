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
    filter(grepl('PD-m|52|NT|RNase',samplename)) %>%
    mutate(prep = ifelse(grepl('SRR',samplename),'ssDNA-seq (ref 1)','TGIRT-seq')) %>%
    filter(distance < 500) %>%
    group_by(prep,distance) %>% 
    summarize(count = n()) %>% 
    ungroup() %>% 
    group_by(prep) %>% 
    do(data_frame(distance = .$distance, 
                  count = .$count/sum(.$count)))  %>%
    tbl_df()

p <- ggplot(data = df, aes(x = distance, y = count, color = prep)) +
    geom_line(size=1.4) +
    labs(x = 'Distance to the Nearest Peak (bp)', y='Fraction of Peaks', color=' ') +
    scale_color_manual(values =c('red','blue')) +
    theme(text = element_text(size=20, face='bold')) +
    theme(legend.position = c(0.7, 0.5))
figurename <- str_c(datapath, 'peak_distance.pdf',sep='/') 
ggsave(p,file = figurename)
message('Plotted: ', figurename)