#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
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
    filter(grepl('^75|umi2',samplename)) %>%
    filter(grepl('nextera|UMI|kh|kq|NEB',samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera-XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('SRR',.$samplename) ~ 'Covaris SRR',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
    filter(!grepl('clustered',samplename)) %>%
    mutate(indel_index = negative_index + positive_index ) %>%
    mutate(number_of_indel = num_D + num_I) %>%
    group_by(prep, samplename, indel_index) %>%
    summarize(number_of_indel = sum(number_of_indel)) %>%
    inner_join(indel_count_df) %>%
    mutate(normalized_indel = number_of_indel / count) %>%
    tbl_df 
    
d <- df %>%
    mutate(normalized_indel = log2(normalized_indel)) %>%
    mutate(indel_index = indel_index) %>%
    tbl_df
    
form <- y ~ poly(x,2)
#form <- y ~ poly(x,1)
source('~/R/legend_to_color.R')
colors <- c('black','salmon','green','orange')
indel_p<-ggplot(data = df %>%
                    filter(grepl('13N|Nextera',prep)) %>%
                    mutate(prep = ifelse(grepl('13N', prep), 'TGIRT-seq',prep)) %>%
                    filter(grepl('umi2|nex', samplename)) %>%
                    mutate(prep = factor(prep, levels = rev(unique(prep)))),
                aes(x = indel_index, y = normalized_indel, color = prep))+
    geom_smooth(se = F,method = 'loess') +
#    geom_smooth(se = F,formula=form, method='lm') +
    geom_point() +
#    geom_line( data= df %>% 
#                   filter(grepl('13N|Nextera',prep)) %>%
#                   mutate(prep = ifelse(grepl('13N', prep), 'TGIRT-seq',prep)) %>%
#                   group_by(prep, indel_index) %>% 
#                   summarize(m = max(normalized_indel)),
#            aes(x=indel_index, color = prep, y = m)) +
    scale_color_manual(values = colors)+
    labs(y = 'Indel freq', color = ' ')+
    scale_x_continuous(breaks = seq(0,10),name='Homopolymer length (nt)') +
    theme(legend.position = c(0.2,0.8)) +
#    theme(legend.text = element_text(size = 25, color = colors)) +
    theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(legend.text = element_text(size=25,face='plain',family = 'Arial')) +
    theme(legend.key.height = unit(2,'line')) 
indel_p <- ggdraw(coloring_legend_text(indel_p)) 
figure_name <- str_c(indel_table_path,'/indel_per_repeat.pdf')
ggsave(indel_p, file = figure_name, width = 7, height = 7)
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