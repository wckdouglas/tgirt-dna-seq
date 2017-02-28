#!/usr/bin/env R

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(purrr)


read_file <- function(filename){
    samplename <- basename(filename)
    samplename <- str_replace(samplename, '.tsv','')
    df <- read_tsv(filename) %>%
        mutate(samplename = samplename)
    return(df)
}

data_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id/bam_files/bed_files'
filenames <- list.files(path = data_path, pattern = '.tsv', full.names = T)

df <- filenames %>% 
    map_df(read_file) %>%
    filter(!grepl('[68]X',samplename)) %>%
    mutate(prep = case_when(grepl('NEB', .$samplename) ~ 'KAPA HiFi + fragmentase',
                            grepl('k[hq]|UMI|7N', .$samplename) ~ 'KAPA HiFi + Covaris',
                            grepl('phus',.$samplename) ~ 'Phusion + Covaris',
                            grepl('q5', .$samplename) ~ 'Q5 + Covaris')) %>%
    tbl_df
