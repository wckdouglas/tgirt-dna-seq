#!/usr/bin/env Rscript

library(stringr)
library(cowplot)
library(readr)
library(dplyr)
library(purrr)


read_file <- function(data_file){
    samplename <- str_split(basename(data_file),'\\.')[[1]][1]
    data_file %>%
    read_csv %>%
    mutate(samplename = samplename) %>%
    return
}

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/bamFiles'
data_files <- list.files(path = datapath, pattern = 'csv', full.names = T)

df <- data_files %>%
    map(read_file)  %>%
    reduce(rbind) %>%
    filter(!grepl('N',base)) %>%
    filter(!is.na(base))

p <- ggplot(data = df, aes(x = positions, y = base_fraction, color = samplename)) +
    geom_line() +
    facet_grid(base~read_end) +
    panel_border()