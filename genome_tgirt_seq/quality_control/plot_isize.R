#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(cowplot)


bampath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/bamFiles/figures'
isize_table <- bampath %>%
	str_c('/insertSize.tsv') %>%
	read_tsv() %>%
	group_by(samplename) %>%
	do(data_frame(
		isize = .$isize,
		counts = .$counts/ sum(.$counts) * 100
	)) %>%
    ungroup() %>%
    dplyr::filter(!grepl('clustered|sim',samplename)) %>%
    dplyr::filter(grepl('nextera|^K12_kh|^K12_kq', samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('SRR',.$samplename) ~ 'Covaris SRR',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
 	tbl_df 


isize_p <- ggplot(data = isize_table, aes(x = isize, y = counts)) +
	geom_line(aes(color = prep, group = samplename), size=1.3) +
	theme(text = element_text(size = 20)) +
	theme(axis.text = element_text(size = 18)) +
	labs(x = 'Insert Size', y = '% of fragments', color = ' ') +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 25, face='bold')) +
    theme(axis.text = element_text(size = 25, face='bold')) +
    xlim(0,600)
    
    
