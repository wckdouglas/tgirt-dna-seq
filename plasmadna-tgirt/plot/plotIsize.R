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
insert_data_path <- stri_c(project_path,'/bamFiles')

rename <- function(x){
    y = ifelse(grepl('PD',x),'TGIRT-ssDNA-seq','ssDNA-seq')
    return(y)
}

insert_df <- insert_data_path %>%
    stri_c('isizeTable.tsv',sep='/') %>%
    read_tsv()%>%
    filter(samplename %in% c('SRR2130052','PD-20X-errorFree'))  %>%
#    filter(chrom == 'Autosomal and Sex Chromosome') %>%
    mutate(samplename = rename(samplename))  %>%
#    group_by(chrom,samplename) %>%
    group_by(samplename) %>%
    do(data_frame(count = .$count/sum(.$count),
                  isize = .$isize)) %>% 
    tbl_df

main_peak <- 167
periodicity <- 10.4
peaks <- main_peak - periodicity * seq(0,6,1)
peaks_df <- data_frame(peak = peaks, colors = c('head',rep('none',length(peaks)-1)))
insert_p <- ggplot(data = insert_df, aes(x=isize, y=count*100)) + 
    geom_bar(stat='identity', color = 'dark blue') + 
    facet_grid(.~samplename)   + 
    theme(axis.text.y = element_text(size=35,face='plain',family = 'Arial')) +
    theme(axis.text.x = element_text(angle=35,hjust=1,size=35, face='plain',family = 'Arial')) +
    theme(text = element_text(size=35, face='bold',family = 'Arial'))+
    labs(x='Insert Size (bp)',y='Percent reads')+
    scale_x_continuous(breaks=seq(0,400,50)) +
    geom_vline(data = peaks_df, aes(xintercept = peak, color = colors), linetype= 2, size=1.5) +
    theme(legend.position='none') +
    scale_color_manual(values = c('green','grey'))