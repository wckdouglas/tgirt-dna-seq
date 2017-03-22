#!/usr/bin/env Rscript

library(stringr)
library(cowplot)
library(extrafont)
library(purrr)
library(dplyr)
library(tidyr)
loadfonts()

setwd('~/tgirt-dna-seq/genome_tgirt_seq/quality_control/')
source('plot_gc.R') #gc_p, #supplemental_p ,lonrenz_curve
source('plot_wgs.R') #wgs_p
source('plot_ends.R') # end_p, prep_end_p
source('plot_errors.R') #error_fig
source('../playing_with_indel/which_base_most_indel.R') #base_indel_p
source('../playing_with_indel/plot_indel.R') #indel_p
source('./plot_sim_end.R') #sim_end_p
source('./plot_sim_isize.R') #sim_size_p
source('./plot_isize.R') #isize_p
source('./plot_bits.R') #small_en_p

figurepath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome//figures'
# main fig plot
p <- ggdraw() +
    draw_plot(wgs_p, 0, 0.5, 0.333, 0.5) +
    draw_plot(gc_p, 0.333, 0.5, 0.333, 0.5) +
    draw_plot(lonrenz_curve, 0.666, 0.5, 0.333, 0.5) +
    draw_plot(error_fig, 0,0,0.5,0.5) +
    draw_plot(indel_p, 0.5,0,0.5,0.5) +
    draw_plot_label(label = letters[1:5], 
                    c(0,0.333,0.666,0, 0.5), c(1,1,1,0.5,0.5), 
                    size = 40, family='Arial') 
figurename <- str_c(figure_path, '/genome.pdf')
ggsave(p, file = figurename, width = 24, height = 18)
message('Plotted:', figurename)
p <- plot_grid(wgs_p, gc_p, lonrenz_curve, 
               error_fig, indel_p, ncol=2,
               labels = letters[1:5],label_size = 40)

figurename <- str_c(figure_path, '/genome.pdf')
ggsave(p, file = figurename, width = 22, height = 20)
message('Plotted:', figurename)




# supplementary genome
p <- ggdraw() +
    draw_plot(isize_p,0,0,0.5,0.5) +
    draw_plot_label(letters[1:3], 
                    x = c(0, 0.5,0), 
                    y = c(1, 1, 0.5), size = 40, family='Arial')
figurename <- str_c(figurepath, '/genome_supplemental_1.pdf')
ggsave(p, file = figurename, width = 15, height = 15)
message('plotted: ', figurename)


# ligation supplement
p<-ggdraw()+
    draw_plot(end_p, 
              0, 0.66, 1, 0.33) +
    draw_plot(prep_end_p, 0, 0.33, 1, 0.33) +
    draw_plot(small_en_p, 0, 0, 0.5, 0.33) +
    draw_plot_label(letters[1:3],c(0,0,0),c(1,0.66,0.33), size =40, family='Arial') 
figurename <- str_c(figurepath, '/genome_supplemental_2.pdf')
ggsave(p, file = figurename, width = 20, height = 20)
message('plotted: ', figurename)


# simulation plots
p <- ggdraw() +
    draw_plot(supplemental_p, 0, 0.5, 0.5,0.5)  +
    draw_plot(sim_size_p, 0,0,0.5,0.5) + 
    draw_plot(sim_end_p, 0.5,0, 0.5, 1) +
    draw_plot_label(letters[1:3], 
                    x = c(0,0, 0.5), 
                    y = c(1,0.5,1),
                    size = 40,
                    family='Arial')
figurename <- str_c(figurepath, '/simulation.pdf')
ggsave(p, file = figurename, width = 25, height = 15)
message('plotted: ', figurename)
    