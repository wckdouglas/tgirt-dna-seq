#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringi)
library(tidyr)
library(cowplot)
library(FBN)
library(extrafont)
library(zoo)
loadfonts()


#source('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/plotDinucleotides.R')
source('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/plotIsize.R')
message('Plotted isize')
source('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/plot_wps_intersample.R')
message('plotted wps intersample')
project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/'
wps_data_path <- stri_c(project_path,'/wpsCTCF')
figure_path <- stri_c(project_path,'/figures')
figurename <- stri_c(figure_path,'/wps_distribution_merged.pdf')

label_rename <- function(x){
    y = ifelse(grepl('^P|^TGIRT',x),'TGIRT-seq','ssDNA-seq')
    return(y)
}

ctcf_df <- wps_data_path %>%
    stri_c('CTCFwps.tsv',sep='/') %>%
    read_tsv()   %>%
    filter(samplename %in% c('SRR2130052','PD_merged'))  %>%
    group_by(samplename, type) %>%
    do(data_frame(
        scale_wps = c(scale(.$wps)),
        position = .$position
    )) %>%
    ungroup() %>%
    mutate(samplename = label_rename(samplename))  %>%
    tbl_df

wps_p <- ggplot(ctcf_df) +
    geom_line(aes(x=position, y = scale_wps, color=samplename), size=1.5) +
    facet_grid(type~., scale='free_y') +
    labs(y= ' ', x ='Position relative to CTCF binding sites', color = ' ')+
    theme(text = element_text(size=35, face='bold', family = 'Arial'))+
    scale_x_continuous(breaks=seq(-1000,1000,200)) +
    theme(axis.text.x = element_text(size=35, angle=50, hjust=1, face='plain',family = 'Arial'))+
    theme(axis.text.y = element_text(size=30,face='plain',family = 'Arial')) +
    theme(legend.position = 'none') +
    scale_color_manual(values = c('black','salmon'))
wps_p <- ggdraw() +
    draw_plot(wps_p, 0.02,0,0.99,1) +
    draw_plot_label('Scaled WPS', 0, 0.37, angle = 90, size = 35)

message('Start plotting')
p <- ggdraw() + 
    draw_plot(insert_p_merge, x= 0.05,y=0.79, width = 0.95, height=0.18) +
    draw_plot(wps_p, x= 0.02,y=0.23, width = 0.98, height=0.5)  +
    #draw_plot(dinucleotide_p + 
    draw_plot(nucleo_p + 
                theme(strip.background = element_blank()) + 
                theme(strip.text.x = element_blank()), x= 0.02,y=0.0, width = 0.98, height=0.21)  +
        draw_plot_label(c('(a)','(b)','(c)'),x = c(0,0,0), y =c(1,0.77,0.30), size=35,family = 'Arial')

ggsave(p, filename = figurename,height = 24,width=14)
message('Plotted:', figurename)
