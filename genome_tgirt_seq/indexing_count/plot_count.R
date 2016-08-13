#!/usr/bin/env Rscript

library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
library(purrr)
library(dplyr)
library(parallel)
library(tidyr)

working_path <- '/stor/work/Lambowitz/cdw2854/jurkatCells/index_tables'
df <- working_path %>%
    str_c('barcode_count.tsv',sep='/') %>%
    read_tsv() %>%
    filter(single_count >4) %>%
    filter(double_count >4) %>%
    tbl_df
    
low_color <- 'yellow'
high_color <- 'red'
p <- ggplot(data = df, aes(x = single_count, y = double_count, fill=log2(count+10))) +
    geom_raster(interpolate = T) +
    labs(x = 'Single barcode family member count', y='Double barcode family member count', fill='log2(count)')  +
    theme(panel.background = element_rect(colour = low_color,fill=low_color)) +
    theme(text = element_text(size=18))+
    scale_fill_gradient(high= high_color,low=low_color) +
    xlim(0,300) +
    ylim(0,300)
figurename <- str_c(working_path,'/index_count.jpg')
ggsave(p, file=figurename,height=5,width=9)
message('plottted ', figurename)
