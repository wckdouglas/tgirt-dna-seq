#!/usr/bin/env Rscript

library(cowplot)
library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)


read_file <- function(tablename){
    samplename <- str_replace(tablename, '.tsv','')
    df <- table_dir %>%
        str_c(tablename, sep='/') %>%
        read_tsv() %>%
        mutate(samplename = samplename) %>%
        filter(strand == '+') %>%
        filter(start > 174) %>%
        filter(start < 370) %>%
        filter(chrom == 'ENST00000621411')  %>%
#        filter(chrom %in% c("ENST00000613854","ENST00000612966","ENST00000621411",
#                            "ENST00000616365","ENST00000359303","ENST00000356476",
#                            "ENST00000614911","ENST00000618052")) %>%
        tbl_df
    return (df)
}

transition <- c('A -> G','G -> A','C -> T','T -> C')

assignTemplate <- function(x){
    ifelse(grepl('^H',x),'RNA','DNA')
}

assignMethod <- function(x){
    ifelse(grepl('bayesian',x),'New clustering','Traditional clustering')
}

table_dir <- '/stor/work/Lambowitz/cdw2854/target-seq/base_tables'
tables <- list.files(path = table_dir, pattern = '.tsv')
merge_df <- tables %>%
	map(read_file) %>%
	reduce(rbind)  %>%
	filter(grepl('.4$',samplename)) %>%
	mutate(template = assignTemplate(samplename)) %>%
	mutate(method = assignMethod(samplename)) %>%
	tbl_df

cov_p <- ggplot(data =merge_df, aes(x=start, y = coverage)) +
    geom_bar(stat='identity',fill='grey') +
    facet_grid(template~method, scale = 'free_y') +
    theme(text = element_text(size = 20, face='bold')) +
    theme(axis.text = element_text(size = 18, face='bold')) +
    labs(x = 'Position on HIST1H3B', y = 'High Quality Base Coverage')
figurename <- str_c(table_dir, '/method_coverage.pdf')
ggsave(cov_p, file = figurename)
message('plotted ', figurename)


