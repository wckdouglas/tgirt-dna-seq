#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(cowplot)
library(tidyr)

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA'
datapath <- str_c(project_path,'/nucleotidesAnaylsis/endsNucleotides')
figurepath <- str_c(project_path,'/figures')

files <- list.files(path = datapath, pattern = '.tsv')
readFile <- function(filename, datapath){
	df <- datapath %>%
		str_c(filename, sep='/') %>%
		read_tsv() %>%
		mutate(filename = filename) %>%
		mutate(name = sapply(filename, function(x) str_split_fixed(x,'_',1))) %>%
		return
}
df <- (files) %>%
    map_df(readFile, datapath) %>%
	select(-total) %>%
	gather(nucleotide, fraction, -index:-name) %>%
	filter(nucleotide!='N') %>%
	tbl_df

p <- ggplot(data = df, aes(x = index, y = fraction, color= nucleotide)) +
	geom_line() +
	facet_grid(End~name)
figurename <- str_c(figurepath, '/nucleotideEnd.pdf')
ggsave(p, file= figurename, width = 15, height = 10)
message('Plotted ', figurename)
