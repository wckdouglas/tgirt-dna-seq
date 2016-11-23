#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)

read_file <- function(filename){
    enzyme <- str_split(basename(filename),'_')[[1]][2]
    filename %>%
        read_tsv(col_names = c('category','value')) %>%
        mutate(enzyme = enzyme) %>%
        filter(category == 'PERCENT_DUPLICATION') %>%
        return
}

rename_enzyme <- function(x){
    if (x == 'kh'){
        'KAPA HIFI'
    }else if (x == 'kq') {
        'KAPA QUANT'
    }else if (x == 'phus'){
        'Phusion'
    }else if (x == 'q5'){
        'Q5'
    }else{
        x
    }
}

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/picard_results'
figure_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures'
figurename <- str_c(figure_path , '/enzyme_duplicates.pdf')
tables <- list.files(path = datapath, pattern = '.duplicates$', full.names=T)
df <- tables %>%
    map(read_file) %>%
    reduce(rbind) %>%
    filter(!grepl('X|S1|Ecoli', enzyme)) %>%
    mutate(enzyme = sapply(enzyme,rename_enzyme)) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(enzyme) %>%
    summarize(
        average_dup = mean(value),
        min_dup = min(value),
        max_dup = max(value)) %>%
    ungroup

p <- ggplot(data = df, aes(x = enzyme, y = average_dup, fill = enzyme)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(ymin = min_dup, ymax = max_dup), width = 0.25) +
    theme(legend.position = 'none') +
    labs(x = ' ', y = 'Average Duplication Rate')
ggsave(p, file = figurename)