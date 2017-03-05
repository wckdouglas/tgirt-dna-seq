#!/usr/bin/env Rscript

library(stringr)
library(cowplot)

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

# main fig plot
ps <- list(wgs_p, gc_p, end_p, error_fig, indel_p, base_indel_p)
ps <- map(ggdraw, ps)
p <- plot_grid(plotlist = ps, 
               ncol=2, 
               labels = letters[1:length(ps)], label_size = 35) 

figurename <- str_c(figure_path, '/genome.pdf')
ggsave(p, file = figurename, width = 22, height = 23)
message('Plotted:', figurename)



# simulation plots
p <- ggdraw() +
    draw_plot(sim_end_p, 0,0, 0.5, 1) +
    draw_plot(sim_size_p, 0.5,0.5,0.5,0.5) + 
    draw_plot(supplemental_p, 0.5, 0, 0.5,0.5)  +
    draw_plot_label(letters[1:3], 
                    x = c(0,0.5,0.5), 
                    y = c(1,1,0.5),
                    size = 35)
figurepath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome//figures'
figurename <- str_c(figurepath, '/simulation.pdf')
ggsave(p, file = figurename, width = 25, height = 15)
message('plotted: ', figurename)

# supplementary genome
p <- ggdraw() +
    draw_plot(isize_p, 0,0,0.333,0.5) +
    draw_plot(prep_end_p, 0.333,0.5,0.667,0.5) +
    draw_plot(small_en_p, 0.333, 0, 0.333,0.5) +
    draw_plot(lonrenz_curve, 0.666,0,0.337,0.5) +
    draw_plot_label(letters[1:5], 
                    x = c(0, 0, 0.333, 0.333, 0.666), 
                    y = c(1, 0.5, 1, 0.5, 0.5),
                    size = 35)
figurepath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome//figures'
figurename <- str_c(figurepath, '/genome_supplemental.pdf')
ggsave(p, file = figurename, width = 28, height = 15)
message('plotted: ', figurename)