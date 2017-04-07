#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(cowplot)
library(FBN)
library(extrafont)
loadfonts()
library(purrr)



source('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/plotDinucleotides.R')
source('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/plotIsize.R')
message('Plotted isize')
source('/stor/home/cdw2854/tgirt-dna-seq/plasmadna-tgirt/plot/plot_wps_intersample.R')
message('plotted wps intersample')
project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/'
wps_data_path <- str_c(project_path,'/wpsCTCF')

read_file <- function(filename){
    read_csv(filename)
}

ctcf_df <- wps_data_path %>%
    list.files(path=., pattern='.csv', full.names=T) %>%
    .[grepl('SRR|P13_mix_UMI', .)]  %>%
    .[grepl('0051|0052|unique', .)]  %>%
    .[grepl('rmdup|unique', .)]  %>%
#    filter(grepl('SRR|P1022', samplename))  %>%
#    .[grepl('0052|umi',.)]  %>%
#    .[grepl('rmdup|umi', .)]  %>%
    map_df(read_csv)   %>%
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
    labs(y= 'Scaled WPS', x ='Position relative to CTCF binding sites (bp)', color = ' ')+
    theme(text = element_text(size=30, family = 'Arial'))+
    scale_x_continuous(breaks=seq(-1000,1000,200)) +
    theme(axis.text.x = element_text(size=30, angle=50, hjust=1, face='plain',family = 'Arial'))+
    theme(axis.text.y = element_text(size=30,face='plain',family = 'Arial')) +
    theme(legend.position = 'none') +
    #theme(legend.key.height = unit(2,'line')) +
    #theme(legend.position = c(0.8,0.3)) +
    scale_color_manual(values = c('salmon','black')) +
    theme(panel.spacing = unit(2, "lines"))
source('~/R/legend_to_color.R')
#wps_p <- ggdraw(coloring_legend_text(wps_p))
wps_p <- ggdraw() +
    draw_plot(wps_p, 0.02,0,0.99,1) #+
    #draw_plot_label('Scaled WPS', 0, 0.37, angle = 90, size = 35)

message('Start plotting')
p <- ggdraw() + 
    draw_plot(insert_p_merge, 0.02,0.42,0.48,0.56) +
    draw_plot(wps_p, 0.5,0.42,0.48,0.56)  +
    draw_plot(dinucleotide_p, 0, 0, 1, 0.38)  +
    draw_plot_label(letters[1:3],
                    x = c(0,0.5,0), 
                    y =c(1,1,0.42), 
                    size=40,family = 'Arial') 

figure_path <- str_c(project_path,'/figures')
figurename <- str_c(figure_path,'/wps_distribution_merged.pdf')
ggsave(p, filename = figurename,height = 16,width=25)
message('Plotted:', figurename)
