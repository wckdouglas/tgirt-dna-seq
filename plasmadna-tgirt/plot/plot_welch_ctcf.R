#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

data_table <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/wpsCTCF/CTCFwps.tsv'
df <- data_table %>%
    read_tsv() %>%
    filter(!grepl('PD3|clustered|subsampled',samplename )) %>%
    filter(grepl('PD_merged|SRR2130051',samplename)) %>%
    mutate(prep = ifelse(grepl('PD',samplename),'TGIRT-seq','ssDNA-seq')) %>%
    group_by(position, prep, type) %>%
    summarize(wps = sum(corrected_wps)) %>%
    ungroup() %>%
    filter(!is.na(wps)) %>%
    group_by(prep,type) %>%
    do(data_frame(
        scale_wps = c(scale(.$wps)),
        position = .$position
    )) %>%
    tbl_df

text_size = 30
p <- ggplot(data = df, 
            aes(color = prep, x = position, y = scale_wps)) +
    geom_line(size = 1.5, alpha=0.7) +
    facet_grid(type~., scale='free_y') +
    labs(x = 'Distance relative to CTCF binding sites',
         y = 'Scaled WPS',
         color = ' ') +
    theme(strip.text.y = element_text(size = text_size, face='bold')) +
    theme(axis.text.x = element_text(size = text_size, face='bold')) +
    theme(axis.text.y = element_text(size = text_size, face='bold')) +
    theme(axis.title = element_text(size = text_size, face = 'bold')) +
    theme(legend.text = element_text(size = text_size, face = 'bold'))+ 
    theme(legend.key.size = unit(.85, "cm"))+ 
    scale_color_manual(values = c('red','black')) +
    theme(legend.position = c(0.8,0.4)) +
    xlim(-900,900) 
figurename <- str_c(dirname(data_table),'/CTCF_binding.pdf')
ggsave(p , file =  figurename, width = 10,height=10)
message('Plotted: ', figurename)