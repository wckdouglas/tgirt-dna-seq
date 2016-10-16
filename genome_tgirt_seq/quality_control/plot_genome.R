#!/usr/bin/env Rscript

library(stringr)
library(cowplot)

source('plot_isize.R')
source('plot_gc.R')
source('plot_wgs.R')

figurename <- str_c(figure_path, '/genome.pdf')
p <- ggdraw() +
	draw_plot(wgs_p, 0.03, 0, 0.47, 0.99) +  
	draw_plot(gc_p, 0.53, 0, 0.47, 0.99) +
	draw_plot_label(label = c('(b)','(c)'), 
	                x= c(0,0.49), y = c(1,1),
	                size = 20)
p <- ggdraw() +
	draw_plot(isize_p, 0.04, 0.5, 0.96, 0.5) +
	draw_plot(p, 0, 0, 1, 0.5) +
	draw_plot_label(label = '(a)',
					x = 0, y = 1, size = 20)
ggsave(p, file = figurename, width = 11, height = 8)
message('Plotted:', figurename)
