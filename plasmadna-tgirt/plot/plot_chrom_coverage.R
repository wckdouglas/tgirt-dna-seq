#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(cowplot)
library(ggrepel)

fai <- '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA' %>%
    str_c('genome_rDNA.fa.fai',sep='/') %>%
    read_tsv(col_names = c('chrom','chrom_length','x','y','z'), 
             col_type = 'cnnnn') %>%
    select(chrom, chrom_length) %>%
    tbl_df

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/bedFiles'
df <- datapath %>%
    str_c('chrom_count.txt', sep='/') %>%
    read_delim(delim=' ', 
               col_names = c('chrom','count','samplename'),
               col_type = 'cnc') %>% 
    inner_join(fai) %>%
    mutate(samplename = str_replace(samplename, '.count','')) %>%
    group_by(samplename) %>%
    do(data_frame(
        count = .$count / sum(.$count),
        chrom = .$chrom,
        chrom_length = .$chrom_length
    )) %>%
    ungroup() %>%
    mutate(method = case_when(grepl('SRR',.$samplename) ~ 'ssDNA-seq',
                              grepl('^P1113|PD[34]', .$samplename) ~ 'TGIRT-Plasma DNA (bad)',
                              grepl('^P1[02]', .$samplename) ~ 'TGIRT-Plasma DNA (Good)',
                              grepl('NT|RN', .$samplename) ~ 'TGIRT-seq Yidan'
                              )
           ) %>%
    tbl_df

model_df <- df %>% 
    group_by(method) %>% 
    nest() %>% 
    mutate(model = map(data, ~lm(count~chrom_length,data=.))) %>% 
    mutate(R = map(model,glance)) %>% 
    unnest(R)  %>%
    select(`r.squared`, method)

plot_df <- df %>% 
    inner_join(model_df) %>%
    mutate(method = str_c(method, '(R-sqrd: ',signif(`r.squared`,3),')'))

p <- ggplot(data = plot_df, 
       aes(color = method, x = chrom_length, y = count)) +
    geom_text(aes(label = samplename)) +
    labs(color = ' ', x = 'Chromosome Length', y = 'Normalized count') +
    geom_smooth(se = F, method='lm') 

model_p <- df %>% 
    group_by(samplename, method) %>% 
    nest() %>% 
    mutate(model = map(data, ~lm(count~chrom_length,data=.))) %>% 
    mutate(R = map(model,glance)) %>% 
    unnest(R)  %>%
    select(`r.squared`, samplename, method) %>%
    set_names(c('r2', 'samplename','method')) %>%
    ggplot(aes(x=method, y =r2, color=method)) +
        geom_jitter() +
        geom_text_repel(aes(label = samplename)) 
    