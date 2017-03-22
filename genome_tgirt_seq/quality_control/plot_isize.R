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
		counts = .$count/ sum(.$count) * 100
	)) %>%
    ungroup() %>%
    dplyr::filter(!grepl('clustered|sim',samplename)) %>%
    dplyr::filter(grepl('nextera|^K12_kh|^K12_kq|UMI', samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('SRR',.$samplename) ~ 'Covaris SRR',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
    filter(grepl('UMI',.$samplename)) %>%
    filter(grepl('_umi2id',.$samplename)) %>%
  	tbl_df 

median_size <- isize_table %>%
    group_by(isize) %>% 
    summarize(counts = sum(counts)) %>% 
    ungroup() %>% 
    filter(counts == max(counts)) %>%
    .$isize

isize_p <- ggplot(data = isize_table, aes(x = isize, y = counts)) +
	geom_line(aes(color = prep, group = samplename), size=1.3) +
	theme(text = element_text(size = 20)) +
	theme(axis.text = element_text(size = 18)) +
	labs(x = 'Insert Size', y = '% fragments', color = ' ') +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 25, face='bold')) +
    theme(axis.text = element_text(size = 25, face='bold')) +
    xlim(0,600)+
    geom_segment(y = 0.78,x = 200, yend = 0.73, xend = 156,
                 arrow = arrow(length = unit(0.5, "cm"))) +
    annotate(geom='text', x = 270, y = 0.78, label = str_c(median_size,' nt'), size = 12, fontface='bold') 

figurepath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures'
figurename <- str_c(figurepath, '/k12_isize.pdf')
ggsave(isize_p, file = figurename)
message('Plotted: ', figurename)
    
    
