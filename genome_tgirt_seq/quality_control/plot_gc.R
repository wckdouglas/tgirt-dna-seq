#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(purrr)
library(stringr)
library(RColorBrewer)

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA'
data_path <- str_c(project_path, 'gc_length',sep='/')
figure_path <- str_c(project_path, 'figures',sep='/')
files <- list.files(path = data_path, pattern = 'DB', full.names = T)

df <- files %>%
    map(read_tsv) %>%
    reduce(rbind) %>%
    mutate(gc = round((gc * 100),digits=0)) %>%
    mutate(seq_lengtg  = round(seq_length,digits=1)) %>%
    group_by(seq_length, gc) %>%
    summarize(count = sum(count)) %>%
    tbl_df

colors <- rev(brewer.pal(9,"Spectral"))
low_color <- colors[1]
p <- ggplot(data = df, aes(x= seq_length, y = gc)) + 
    geom_raster(aes(fill = log2(count)),interpolate = T) +
    scale_fill_gradientn(colors = colors) +
	labs(y = 'GC %',x = 'Insert size', color = 'log2(Count)')+
	theme(panel.background = element_rect(colour = low_color,fill=low_color))
figurename <- str_c(figure_path, '/gc_length.jpg')
ggsave(p, file = figurename, height=5,width=7)
message('Plotted ' , figurename)
