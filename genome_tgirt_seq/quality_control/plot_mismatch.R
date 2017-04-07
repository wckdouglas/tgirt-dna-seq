#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(cowplot)


read_file <- function(tablename){
    samplename <- str_replace(basename(tablename),'.xls','')
    df <- read_tsv(tablename) %>%
        gather(mismatch, count, -read_pos, -sum) %>%
        separate(col = mismatch, into = c('from','to'),sep = '2') %>%
        group_by(from, to) %>%
        summarize(count = sum(count)) %>%
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
    mutate(mismatch_rate = count/10000) %>%
    mutate(mismatch = str_c(from, to, sep='-to-')) %>%
    filter(grepl('75bp|clustered_family|umi2id.M',samplename)) %>%
    mutate(prep = case_when(
                    grepl('nextera', .$samplename) ~ 'Nextera-XT',
                    grepl('clustered_family',. $samplename) ~ 'TGIRT-seq\n(Error-corrected)',
                    TRUE ~ 'TGIRT-seq')) %>%
    group_by(mismatch, prep, from) %>%
    summarize(error = mean(mismatch_rate),
              deviation = sd(mismatch_rate)) %>%
    tbl_df


mismatch_p <- ggplot(data = df, aes(x = mismatch, y = error, fill = prep)) +
    geom_bar(stat='identity', position = 'dodge') +
    geom_errorbar(position = position_dodge(width = 1), width=0.2,
                  aes(ymin = error - deviation, ymax = error+deviation)) +
    facet_grid(prep~from, scale='free_x') +
    labs(fill = ' ', x= ' ', y = 'Error rate per\nerror containing read') +
    theme(axis.text.x = element_text(angle = 60, size=30, family='Arial', hjust = .5, vjust = .5)) +
    theme(axis.text.y = element_text(size=30, family='Arial')) +
    theme(text = element_text(size=30, family='Arial')) +
    theme(legend.key.height = unit(2,'line')) +
    theme(legend.position = 'none')+
    scale_fill_manual(values = c('salmon','black','green4')) +
    panel_border()
figurename <- str_c(table_path, '/genome_mismatch.pdf')
ggsave(mismatch_p, file = figurename)
message('plotted: ', figurename)