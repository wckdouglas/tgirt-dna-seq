#!/usr/bin/env Rscript

library(cowplot)
library(stringr)

setwd('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/')
source('plot_wps_intersample.R') #dist_p
source('plot_wps_peak_distance.R') #nucleo_p
source('plot_tissue_cor.R') #nucleo_p

p <- ggdraw()+
    draw_plot(dist_p, 0.01,0.5,0.49,0.48) + 
    draw_plot(nucleo_p, 0.52,0.46,0.48,0.52) + 
    draw_plot(tissue_cor_p, 0.05,0,0.9,0.48) + 
    draw_plot_label(letters[1:3], 
                    x = c(0,0.5,0), 
                    y = c(1,1,0.5),
                    size=50,family = 'Arial') 
figurepath <-  '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures'
figurename <- str_c(figurepath, '/nucleosome.pdf')
ggsave(p, file = figurename, height=15, width=18)
message(figurename)
