#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(cowplot)

bampath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/bamFiles'
isize_table <- bampath %>%
	str_c('/isizeTable.tsv') %>%
	read_tsv() %>%
	group_by(samplename) %>%
	do(data_frame(
		isize = .$isize,
		count = .$count/ sum(.$count) * 100
	)) %>%
    ungroup() %>%
    dplyr::filter(!grepl('clustered|sim',samplename)) %>%
	tbl_df 


isize_p <- ggplot(data = isize_table, aes(x = isize, y = count)) +
	geom_line(aes(color = samplename), size=1.3) +
	theme(text = element_text(size = 20)) +
	theme(axis.text = element_text(size = 18)) +
	labs(x = 'Insert Size', y = 'Percentage of fragments', color = ' ')
    
