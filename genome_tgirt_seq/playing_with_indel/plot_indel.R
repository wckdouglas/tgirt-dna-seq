#!/usr/bin/env Rscript

library(readr)
library(purrr)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(cowplot)

indel_index_table <- '/stor/work/Lambowitz/ref/Ecoli/k12_mg1655_repeat_index.tsv'
indel_index_df <- read_tsv(indel_index_table) %>%
    mutate(indel_index=negative_index + positive_index) %>% 
    tbl_df
indel_count_df <-  indel_index_df %>% 
    group_by(indel_index) %>% 
    summarize(count = n())

read_indel_table <- function(tablename){
    samplename <- basename(str_replace(tablename,'.tsv',''))
    tablename %>%
    read_tsv() %>%
    mutate(samplename = samplename) %>%
    tbl_df %>%
    return()
}


indel_table_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/indel_table'
df <- list.files(path = indel_table_path, pattern = '.tsv', full.names = T) %>%
    map_df(read_indel_table) %>%
    inner_join(indel_index_df)%>% 
    select(grep('samplename|num_D|num_I|index',names(.))) %>%
    filter(!grepl('MiSeq|Ecoli',samplename)) %>%
    mutate(prep = ifelse(grepl('nextera',samplename),'Nextera XT',
                         ifelse(grepl('clustered',samplename),'Clustered TGIRT-seq','TGIRT-seq'))) %>%
    mutate(indel_index = negative_index + positive_index ) %>%
    mutate(number_of_indel = num_D + num_I) %>%
    group_by(prep, samplename, indel_index) %>%
    summarize(number_of_indel = sum(number_of_indel)) %>%
    inner_join(indel_count_df) %>%
    mutate(normalized_indel = number_of_indel / count) %>%
    tbl_df
    
form <- y ~ poly(x,1, raw=T)
p<-ggplot(data = df, aes(x = indel_index, y = normalized_indel, color = prep))+
    geom_smooth(method='lm', formula =form, se = F) +
    geom_point() +
    ggpmisc::stat_poly_eq(aes(label = ..eq.label..), 
                 formula = form, 
                 label.x.npc = "center",
                 label.y.npc = 'top',
                 parse = TRUE) +
    ggpmisc::stat_poly_eq(aes(label = ..rr.label..),
                 formula = form, 
                 label.y.npc = "center",
                 label.x.npc = 'center',
                 parse = TRUE)+
    labs(x = 'Repeat Index on Genome', y = 'Average Indel per Read', color = ' ')+
    scale_x_continuous(breaks = seq(0,10),name='Repeat units') +
    theme(legend.position = c(0.2,0.8))
figure_name <- str_c(indel_table_path,'/indel_per_repeat.pdf')
ggsave(p, file = figurename)
message('Plotted: ', figurename)