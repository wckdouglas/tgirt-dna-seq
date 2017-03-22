#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringi)
library(tidyr)
library(cowplot)
library(FBN)
library(extrafont)
loadfonts()

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/'
insert_data_path <- stri_c(project_path,'/insertSize')

rename <- function(x){
    y = ifelse(grepl('^P',x),'TGIRT-seq','ssDNA-seq')
    return(y)
}

insert_df <- insert_data_path %>%
    stri_c('isizeTable.tsv',sep='/') %>%
    read_tsv()%>%
    filter(grepl('SRR2130051|SRR2130052|P1022',samplename))  %>%
    filter(grepl('umi2id|SRR', samplename)) %>%
    filter(!grepl('cluster',samplename))  %>%
    filter(isize > 22) %>%
    mutate(prep = rename(samplename))  %>%
    group_by(isize,samplename, prep) %>%
    summarize(count = sum(count)) %>%
    ungroup() %>%
    group_by(prep,samplename) %>%
    do(data_frame(count = .$count/sum(.$count),
                  isize = .$isize)) %>% 
    ungroup() %>%
    mutate(prep = factor(prep, level = rev(unique(prep)))) %>%
    arrange(isize) %>%
    tbl_df

main_peak <- insert_df %>% filter(count == max(count)) %>% .$isize
periodicity <- 10.4
peaks <- main_peak - periodicity * seq(0,10,1)
peaks_df <- data_frame(peak = peaks, colors = c('head',rep('none',length(peaks)-1)))


insert_p_merge <- ggplot(data = insert_df) + 
    geom_line(size=2,alpha=0.8,aes(x=isize, y=count*100, 
                                        color = prep,
                                        group = samplename)) +
    theme(axis.text.y = element_text(size=35,face='plain',family = 'Arial')) +
    theme(axis.text.x = element_text(angle=50,hjust=1,size=30, face='plain',family = 'Arial')) +
    theme(text = element_text(size=35, family = 'Arial'))+
    labs(x='Fragment length (nt)',y='Percent reads', color = ' ')+
    theme(legend.key.size = unit(11,'mm')) +
    theme(legend.position = c(0.8,0.9)) +
    scale_x_continuous(breaks=seq(0,401,50), limits=c(0,400)) +
    geom_vline(data = peaks_df, aes(xintercept = peak), 
               linetype= 2, size=1, color = 'gray1', alpha = 0.5) +
    theme(legend.position = c(0.8,0.4)) +
    geom_segment(y = 2.4,x = 200, yend = 2.1, xend = 170,
                 arrow = arrow(length = unit(0.5, "cm"))) +
    annotate(geom='text', x = 234, y = 2.5, label = '167 nt', size = 12, fontface='bold') +
    ylim(0,2.5) +
    scale_color_manual(values = c('black','salmon'))
source('~/R/legend_to_color.R')
insert_p_merge <- ggdraw(coloring_legend_text(insert_p_merge))
figurename <- str_c(insert_data_path, '/plasma_insert_profile.pdf')
ggsave(insert_p_merge, file = figurename, height = 8)
message('Plotted: ', figurename)
