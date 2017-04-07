#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(data.table)
library(tidyr)
library(stringr)

df <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/mismatch_profiles/mismatch_profiles.tsv' %>%
    read_tsv() %>%
    separate(col = mismatch,into = c('ref','read'), sep='>') %>%
    group_by(samplename, ref) %>%
    do(data_frame(
        read =.$read,
        counts = .$counts,
        fraction = .$counts/sum(.$counts)
    )) %>%
    ungroup() %>%
    filter(ref != read) %>%
    mutate(prep = case_when(
        grepl('nextera',.$samplename) ~ 'Nextera-XT',
        grepl('family', .$samplename) ~ 'TGIRT-seq\n(Error-corrected)',
        grepl('umi2id', .$samplename) ~ 'TGIRT-seq'
    )) %>%
    filter(!is.na(prep)) %>%
    mutate(fraction = fraction * 1e5)
    tbl_df

sub_p<-ggplot(data= df, aes(x = str_c(ref, '-to-', read), y = fraction, color = prep)) +
    geom_jitter() +
    facet_grid(prep~ref, scales = 'free_x') +
    labs(x = 'Substitution', y = expression(Substitution~rate~(x10^5)), color = ' ', parse=T) +
    theme(text = element_text(size=30,  family = 'Arial')) +
    theme(strip.text = element_text(size = 30,  family = 'Arial'))+
    theme(axis.text = element_text(size = 30, family = 'Arial'))+
    theme(axis.text.x = element_text(angle = 70, hjust=0.5, vjust=0.5, family='Arial', size=30)) +
    theme(legend.position='none')+
    panel_border()
