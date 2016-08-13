#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)

projectpath <- '/scratch/cdw2854/plasmaDNA'
datapath <- str_c(projectpath,'/rnaCalls')
figurepath <- str_c(projectpath, '/figures')
figurename <- str_c(figurepath, '/rna_deconvolved.pdf')
p <- datapath%>% 
	str_c('rnaCount.tsv',sep='/') %>%
	read_tsv() %>%
	group_by(samplename) %>%
	do(data_frame(rna = .$rnaType,
				  scores = .$scores/sum(.$scores) * 100)) %>%
	ggplot(aes(x=samplename, y= scores, fill=rna)) +
		geom_bar(stat='identity') +
		theme(axis.text.x = element_text(angle=90,hjust=1, vjust=0.5)) +
		labs(x=' ', y = 'Percentage',fill=' ')
ggsave(p, file= figurename, width=15,height=10)
message('Plotted ', figurename)
