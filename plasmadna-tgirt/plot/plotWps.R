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
figurename <- stri_c(figure_path,'/wps_distribution.pdf')

label_rename <- function(x){
    y = ifelse(grepl('^P|^TGIRT',x),'TGIRT-seq','ssDNA-seq')
    return(y)
}

ctcf_df <- wps_data_path %>%
    stri_c('CTCFwps.tsv',sep='/') %>%
    read_tsv()  


sim_data <- ctcf_df %>%
    filter(grepl('sim',samplename)) %>%
    filter(samplename != 'sim_bed') %>%
    mutate(samplename = label_rename(samplename))%>% 
    dplyr::rename(sim_wps = corrected_wps) %>%
    select(-wps)


ctcf_df <- ctcf_df %>%
    filter(samplename %in% c('SRR2130052','PD_merged'))  %>%
    group_by(samplename) %>%
    do(data_frame(wps = .$wps,
                  corrected_wps = .$corrected_wps,
                  position = .$position,
                  type = .$type)) %>%
    ungroup() %>%
    mutate(samplename = label_rename(samplename)) %>%
    inner_join(sim_data, by=c('samplename','position','type')) %>% 
    mutate(adjusted_wps = corrected_wps - sim_wps)  %>%
    mutate(scaled_wps = scale(adjusted_wps)) 

wpsPlot <- function(sample){
    p <- ctcf_df %>%
        filter(samplename == sample) %>%
        ggplot() +
            geom_line(aes(x=position, y = adjusted_wps), 
                    color='dark blue', size=1.5) +
            geom_line(aes(x = position, y=sim_wps),
                      alpha=0.5, color = 'grey') +
            facet_grid(type~., scale='free_y') +
            labs(x= ' ', y =' ')+
            theme(text = element_text(size=35, face='bold', family = 'Arial'))+
            scale_x_continuous(breaks=seq(-1000,1000,200)) +
            theme(axis.text.x = element_text(size=35, angle=50, hjust=1, face='plain',family = 'Arial'))+
            theme(axis.text.y = element_text(size=35,face='plain',family = 'Arial'))  
    return (p)
}


samples = rev(c('TGIRT-seq','ssDNA-seq'))
wps_ps <- lapply(samples,wpsPlot)
wps_p <- plot_grid(wps_ps[1][[1]], 
                   wps_ps[2][[1]] + theme(axis.text.y = element_blank())) +
    draw_label( 'Distance to CTCF start site (bp)', x = 0.5, y = 0.05, fontface = 'bold',fontfamily='Arial',size = 35)+
    draw_label( 'Scaled WPS', y = 0.6, x = 0.02, fontface = 'bold',fontfamily='Arial',size = 35, angle=90)

message('Start plotting')
p <- ggdraw() + 
    draw_plot(insert_p, x= 0.05,y=0.8, width = 0.95, height=0.199) +
    draw_plot(wps_p, x= 0.02,y=0.3, width = 0.98, height=0.49)  +
    #draw_plot(dinucleotide_p + 
    draw_plot(nucleo_p + 
                theme(strip.background = element_blank()) + 
                theme(strip.text.x = element_blank()), x= 0.02,y=0.0, width = 0.98, height=0.28)  +
        draw_plot_label(c('(a)','(b)','(c)'),x = c(0,0,0), y =c(1,0.8,0.315), size=35,family = 'Arial')

ggsave(p, filename = figurename,height = 20,width=18)
message('Plotted:', figurename)
