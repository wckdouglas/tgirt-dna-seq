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
ps <- list(wgs_p, gc_p, error_fig, indel_p)
p <- plot_grid(plotlist = ps, 
               ncol=2, 
               labels = letters[1:length(ps)], label_size = 35) 

figurename <- str_c(figure_path, '/genome.pdf')
ggsave(p, file = figurename, width = 18, height = 18)
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
    draw_plot(isize_p,0.5,0,0.5,1) +
    draw_plot_label(letters[1:2], 
                    x = c(0, 0.5), 
                    y = c(1, 1), size = 35)
figurepath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome//figures'
figurename <- str_c(figurepath, '/genome_supplemental_1.pdf')
ggsave(p, file = figurename, width = 15, height = 7)
message('plotted: ', figurename)


# ligation supplement
ps <- list(small_en_p,
            lonrenz_curve,
            base_indel_p)
p <- plot_grid(plotlist = ps, ncol = 3,
               labels =letters[2:4], 
               label_size = 35)
p<-ggdraw()+
    draw_plot(prep_end_p, 0, 0.5, 1, 0.5) +
    draw_plot(p, 0, 0, 1, 0.5) +
    draw_plot_label(letters[1],0,1, size =35) 
figurepath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome//figures'
figurename <- str_c(figurepath, '/genome_supplemental_2.pdf')
ggsave(p, file = figurename, width = 28, height = 15)
message('plotted: ', figurename)

    