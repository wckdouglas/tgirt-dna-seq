#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringi)
library(cowplot)
library(tidyr)

datapath <- '/Users/wckdouglas/plasmaDNA/results/plasmaResults/nucleotidesAnalysis'
figurepath <- '/Users/wckdouglas/plasmaDNA/figures'

files <- list.files(path = datapath, pattern = '.tsv')
readFile <- function(filename, datapath){
	df <- datapath %>%
		stri_c(filename, sep='/') %>%
		read_tsv() %>%
		mutate(filename = filename) %>%
		separate(filename, 'name', '_', remove=T) %>%
		return
}
df <- lapply(files, readFile, datapath) %>%
	do.call(rbind,.) %>%
	select(-total) %>%
	gather(nucleotide, fraction, -index:-name) %>%
	filter(nucleotide!='N') %>%
	tbl_df

p <- ggplot(data = df, aes(x = index, y = fraction, color= nucleotide)) +
	geom_line() +
	facet_grid(End~name)
figurename <- stri_c(figurepath, '/nucleotideEnd.pdf')
ggsave(p, file= figurename, width = 15, height = 10)
message('Plotted ', figurename)
