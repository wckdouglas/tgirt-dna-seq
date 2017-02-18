#!/usr/bin/evn Rscript

library(readr)
library(dplyr)
library(tibble)
library(purrr)

read_table <- function(filename) {
    read_tsv(filename,skip = 6) %>% 
    head(1) %>% 
    gather() %>%
    set_names(c('metric','value')) %>%
    mutate(samplename = str_split(basename(filename),'\\.')[[1]][1]) %>%
    mutate(value = as.numeric(value)) %>%
    return()
}

project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
WGS_files <- list.files(path = picard_path, 
                        pattern = '.wgs_metric',
                        full.names = T)
df <- WGS_files %>%
    map(read_table) %>%
    reduce(rbind)%>%
    filter(grepl('^PCT_[0-9]+X$',metric)) %>%
    filter(grepl('K12_k[qh]|nextera',samplename)) %>%
    filter(!grepl('sim',samplename)) %>%
    mutate(prep = ifelse(grepl('nextera',samplename),'Nextera XT',
                         ifelse(grepl('cluster',samplename),'Clustered TGIRT-seq','TGIRT-seq')))  %>%
    mutate(coverage = str_replace_all(metric,'PCT_','')) %>%
    mutate(coverage = str_replace_all(coverage,'X','')) %>%
    mutate(coverage = as.numeric(coverage)) %>%
    tbl_df

p <- ggplot(data = df , aes(x = coverage, y = value, 
                            group=samplename, 
                            color= prep)) +
    geom_line()+
    xlim(1,30) + 
    theme(legend.position = c(0.7,0.5)) +
    labs(color = ' ', x = 'Coverage', y = 'Percentage of bases\nat this Coverage or Better')