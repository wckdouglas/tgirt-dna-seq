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
    filter(grepl('SRR2130052|P1',samplename))  %>%
    filter(!grepl('cluster',samplename))  %>%
    filter(isize > 22) %>%
    mutate(samplename = rename(samplename))  %>%
    group_by(isize,samplename) %>%
    summarize(count = sum(count)) %>%
    ungroup() %>%
    group_by(samplename) %>%
    do(data_frame(count = .$count/sum(.$count),
                  isize = .$isize)) %>% 
    tbl_df

main_peak <- 167
periodicity <- 10.4
peaks <- main_peak - periodicity * seq(0,10,1)
peaks_df <- data_frame(peak = peaks, colors = c('head',rep('none',length(peaks)-1)))
insert_p <- ggplot(data = insert_df, aes(x=isize, y=count*100)) + 
    geom_bar(stat='identity', color = 'dark blue') + 
    facet_grid(.~samplename)   + 
    theme(axis.text.y = element_text(size=35,face='plain',family = 'Arial')) +
    theme(axis.text.x = element_text(angle=50,hjust=1,size=30, face='plain',family = 'Arial')) +
    theme(text = element_text(size=35, face='bold',family = 'Arial'))+
    labs(x='Fragment length (nt)',y='Percent Reads')+
    scale_x_continuous(breaks=seq(0,401,50), limits=c(0,400)) +
    geom_vline(data = peaks_df, aes(xintercept = peak, color = colors), linetype= 2, size=1) +
    theme(legend.position='none') +
    scale_color_manual(values = c('green','grey'))

insert_p_merge <- ggplot(data = insert_df) + 
    geom_line(size = 1.5, alpha=0.8,aes(x=isize, y=count*100, color = samplename)) +
    theme(axis.text.y = element_text(size=35,face='plain',family = 'Arial')) +
    theme(axis.text.x = element_text(angle=50,hjust=1,size=30, face='plain',family = 'Arial')) +
    theme(text = element_text(size=35, face='bold',family = 'Arial'))+
    labs(x='Fragment length (nt)',y='Percent Reads', color = ' ')+
    theme(legend.key.size = unit(11,'mm')) +
    scale_x_continuous(breaks=seq(0,401,50), limits=c(0,400)) +
    geom_vline(data = peaks_df, aes(xintercept = peak), 
               linetype= 2, size=1, color = 'gray1', alpha = 0.5) +
    theme(legend.position = c(0.8,0.5)) +
    geom_segment(y = 2.4,x = 200, yend = 2.1, xend = 170,
                 arrow = arrow(length = unit(0.5, "cm"))) +
    annotate(geom='text', x = 234, y = 2.5, label = '167 nt', size = 12, fontface='bold') +
    ylim(0,2.5) +
    scale_color_manual(values = c('salmon','black'))
figurename <- str_c(insert_data_path, '/plasma_insert_profile.pdf')
ggsave(insert_p_merge, file = figurename, height = 8)
message('Plotted: ', figurename)
