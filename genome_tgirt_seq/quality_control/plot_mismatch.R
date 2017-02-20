#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(cowplot)


read_file <- function(tablename){
    samplename <- str_replace(basename(tablename),'.xls','')
    df <- read_tsv(tablename, skip = 1) %>%
        gather(mismatch, count, -read_pos, -sum) %>%
        separate(col = mismatch, into = c('from','to'),sep = '2') %>%
        group_by(from, to) %>%
        summarize(sum = sum(sum),
                  count = sum(count)) %>%
        mutate(mismatch_rate = count) %>%
        mutate(samplename = samplename) %>%
        ungroup() %>%
        tbl_df
    return (df)
}

table_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/mismatch_profiles'
tablenames <- list.files(path = table_path, 
                         pattern = '.xls', 
                         full.names = T)
df <- tablenames %>%
    map_df(read_file) %>%
    mutate(mismatch = str_c(from, to, sep='-to-')) %>%
    mutate(prep = ifelse(grepl('K12',samplename),'TGIRT-seq','Nextera-XT')) %>%
    mutate(cluster = ifelse(grepl('cluster',samplename),' error-corrected','')) %>%
    mutate(prep = str_c(prep, cluster)) %>%
    group_by(mismatch, prep, from) %>%
    summarize(error = mean(mismatch_rate),
              deviation = sd(mismatch_rate)) %>%
    tbl_df

p <- ggplot(data = df, aes(x = mismatch, y = error, fill = prep)) +
    geom_bar(stat='identity', position = 'dodge') +
    geom_errorbar(position = position_dodge(width = 1), width=0.2,
                  aes(ymin = error - deviation, ymax = error+deviation)) +
    facet_grid(.~from, scale='free_x') +
    labs(fill = ' ', x= ' ', y = 'Error base count') +
    theme(axis.text.x = element_text(angle = 60, size=25, face='bold', hjust = .5, vjust = .5)) +
    theme(text = element_text(size=25, face='bold')) +
    theme(legend.key.height = unit(2,'line')) +
    scale_fill_manual(values = c('light sky blue','salmon','green4'))
figurename <- str_c(table_path, '/genome_mismatch.pdf')
ggsave(p, file = figurename)
message('plotted: ', figurename)