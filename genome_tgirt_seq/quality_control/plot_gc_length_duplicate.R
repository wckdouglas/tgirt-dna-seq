#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(RColorBrewer)

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/gc_length_duplicate'
df <- datapath %>%
    str_c('DB-20X-2-errorFree.tsv',sep='/') %>%
    read_tsv() 

colors <- rev(brewer.pal(9,"Spectral"))
low_color <- colors[1]
p <- ggplot(data =df, aes(x = gc_per, y = (count))) +
    stat_density2d(aes(fill = ..density..),contour = FALSE, geom="tile") +
    scale_fill_gradientn(colors = colors) +
    labs(y = 'Duplicate',x = 'GC %')+
    theme(panel.background = element_rect(fill=low_color))