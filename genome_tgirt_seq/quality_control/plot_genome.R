#!/usr/bin/env Rscript

library(stringr)
library(cowplot)

setwd('~/tgirt-dna-seq/genome_tgirt_seq/quality_control/')
source('plot_isize.R')
source('plot_gc.R')
source('plot_wgs.R')

fix_color<-function(p){
    p + scale_color_manual(values = c('light sky blue','salmon'))
}

figurename <- str_c(figure_path, '/genome.pdf')

p <- ggdraw() +
    draw_plot(fix_color(isize_p), 0.03,0.5,0.46,0.5)+
    draw_plot(fix_color(gc_p), 0.52, 0.5,0.48,0.5) +
    draw_plot(wgs_p, 0, 0, 1, 0.5) +
    draw_plot_label(label = c('(a)','(b)','(c)'), 
                x= c(0, 0, 0.49), y = c(1,0.54, 1),
                size = 20)     

ggsave(p, file = figurename, width = 12, height = 10)
message('Plotted:', figurename)
