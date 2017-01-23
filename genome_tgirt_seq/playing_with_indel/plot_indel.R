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
    filter(!grepl('MiSeq|Ecoli|phus|q5|NEB',samplename)) %>%
    mutate(prep = case_when(grepl('nextera', .$samplename)~'Nextera XT',
                            grepl('clustered', .$samplename)~'Clustered TGIRT-seq',
                            grepl('NEB', .$samplename)~'NEB TGIRT-seq',
                            grepl('K12', .$samplename)~'TGIRT-seq')) %>%
    filter(!grepl('clustered',samplename)) %>%
    mutate(indel_index = negative_index + positive_index ) %>%
    mutate(number_of_indel = num_D + num_I) %>%
    group_by(prep, samplename, indel_index) %>%
    summarize(number_of_indel = sum(number_of_indel)) %>%
    inner_join(indel_count_df) %>%
    mutate(normalized_indel = number_of_indel / count) %>%
    tbl_df
    
form <- y ~ poly(x,2)
p<-ggplot(data = df, aes(x = indel_index, y = normalized_indel, color = prep))+
#    geom_smooth(se = F,method = 'loess') +
geom_smooth(se = F,formula=form, method='lm') +
    geom_point() +
    scale_color_manual(values = c('lightskyblue','salmon'))+
    labs(y = 'Average Indel per Read\nper Mononucleotide Run', color = ' ')+
    scale_x_continuous(breaks = seq(0,10),name='Run Length (nt)') +
    theme(legend.position = c(0.2,0.8)) +
    theme(text = element_text(size = 25, face='bold')) +
    theme(axis.text = element_text(size = 25, face='bold')) +
    theme(legend.key.height = unit(2,'line'))
figure_name <- str_c(indel_table_path,'/indel_per_repeat.pdf')
ggsave(p, file = figure_name, width = 7, height = 7)
message('Plotted: ', figure_name)


library(broom)
dodge <- position_dodge(width=0.9)
model_p <- df %>% 
    group_by(prep)%>% 
    nest() %>% 
    mutate(model = map(data, ~lm(normalized_indel~poly(indel_index,2),data = .))) %>% 
    mutate(model_res = map(model, tidy)) %>%
    unnest(model_res) %>%
    tbl_df %>%
    ggplot(aes(x=term, y = estimate, fill=prep))+
        geom_bar(position=dodge, stat='identity') +
        geom_errorbar(position=dodge, width=0.5,
                      aes(ymax=estimate+std.error, ymin = estimate-std.error))