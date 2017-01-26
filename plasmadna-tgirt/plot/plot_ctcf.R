#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
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
wps_data_path <- str_c(project_path,'/wpsCTCF')
figure_path <- str_c(project_path,'/figures')
figurename <- str_c(figure_path,'/wps_distribution_merged.pdf')


ctcf_df <- wps_data_path %>%
    stri_c('CTCFwps.tsv',sep='/') %>%
    read_tsv()   %>%
    filter(grepl('SRR|P1022', samplename))  %>%
    filter(grepl('0052|umi', samplename))  %>%
    filter(grepl('rmdup|umi', samplename))  %>%
    mutate(prep = case_when(grepl('^P|^TGIRT',.$samplename) ~ 'TGIRT-seq',
                                grepl('SRR',.$samplename) ~ 'ssDNA-seq (ref.2)'))  %>%
    group_by(prep, type, position) %>%
    summarize(wps = sum(wps)) %>%
    ungroup() %>%
    group_by(prep, type) %>%
    do(data_frame(
        scale_wps = c(scale(.$wps)),
        position = .$position,
        wps = .$wps
    )) %>%
    ungroup() %>%
    mutate(type = str_replace(type,'bp','nt')) %>%
    tbl_df

wps_p <- ggplot(ctcf_df) +
    geom_line(aes(x=position, y = scale_wps, color=prep), 
              size=1.5, alpha = 0.6) +
    facet_grid(type~., scale='free_y') +
    labs(y= ' ', x ='Position relative to CTCF binding sites', color = ' ')+
    theme(text = element_text(size=35, face='bold', family = 'Arial'))+
    scale_x_continuous(breaks=seq(-1000,1000,200)) +
    theme(axis.text.x = element_text(size=35, angle=50, hjust=1, face='plain',family = 'Arial'))+
    theme(axis.text.y = element_text(size=30,face='plain',family = 'Arial')) +
#    theme(legend.position = 'none') +
    theme(legend.position = c(0.8,0.5)) +
    scale_color_manual(values = c('salmon','black'))
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
        draw_plot_label(c('(a)','(b)','(c)'),x = c(0,0,0), y =c(1,0.77,0.24), size=35,family = 'Arial')

ggsave(p, filename = figurename,height = 24,width=14)
message('Plotted:', figurename)
