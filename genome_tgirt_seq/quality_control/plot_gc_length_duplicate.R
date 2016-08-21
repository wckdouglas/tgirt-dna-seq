#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(RColorBrewer)
library(purrr)

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/gc_length_duplicate'
files <- list.files(path = datapath, pattern='DB',full.names = T)


df <- files[grep('errorFree',files,invert=T)] %>%
    map(read_tsv) %>%
    reduce(rbind) 

gc_df <- df %>%
    group_by(gc_per) %>%
    summarize(mean_count = mean(count)) 

len_df <- df %>%
    group_by(length) %>%
    summarize(mean_count = mean(count))  

p_gc <- ggplot(data = gc_df, aes(x = gc_per, y = mean_count)) +
    geom_point() +
    labs(x = '%GC content',y = 'Mean duplication number')

p_len <- ggplot(data = len_df, aes(x=length, y = mean_count)) +
    geom_point() +
    labs(x = 'Length', y = ' ')

p <- plot_grid(p_gc, p_len, labels = c('(a)','(b)'))