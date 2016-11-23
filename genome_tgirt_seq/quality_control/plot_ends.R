#!usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)

datapath <- '/stor/scratch/Lambowitz/cdw2854/dna'
files <- list.files(path = datapath, pattern = '.csv')

rename_sample <- function(x){
    ifelse(grepl('NA',x),'NA12878 TGIRT',
           ifelse(grepl('EG',x),'B-strain E.coli TGIRT',x)
           )
}

get_sample_number <- function(x){
    ifelse(grepl('bt|EG2',x),'2','1')
}

read_files <- function(filename){
    sample_name = str_split(filename,'\\.')[[1]][1]
    d <- datapath %>%
        str_c(filename,sep='/') %>%
        read_csv() %>%
        mutate(filename = sample_name)
    return(d)
}

df <- files %>%
    map(read_files) %>%
    reduce(rbind) %>%
    filter(!grepl('clustered',filename)) %>%
    mutate(samplename = rename_sample(filename)) %>%
    mutate(sample_num = get_sample_number(filename)) %>%
    mutate(filename = str_c(samplename, sample_num, sep=' ')) %>%
    mutate(read_end = factor(read_end, levels=c("5'","3'")))

p <- ggplot(data = df, aes(x = positions, color = filename, 
                           y = base_fraction)) +
    geom_line() +
    facet_grid(base~read_end) +
    labs(x = 'Positions',y='Fraction of base',color=' ') +
    panel_border()