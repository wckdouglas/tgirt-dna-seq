#!/usr/bin/env Rscript

library(stringr)
library(cowplot)

source('plot_gc.R')
source('plot_wgs.R')

figurename <- str_c(figure_path, '/genome.pdf')
p <- ggdraw() +
	draw_plot(wgs_p, 0.03, 0, 0.47, 0.99) +  
	draw_plot(gc_p + theme(legend.position = 'none'), 0.53, 0, 0.49, 0.99) +
	draw_plot_label(label = c('(a)','(b)'), 
	                x= c(0,0.49), y = c(1,1),
	                size = 20)
ggsave(p, file = figurename, width = 11)
message('Plotted:', figurename)
