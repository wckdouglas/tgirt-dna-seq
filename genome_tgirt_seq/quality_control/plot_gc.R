#!/usr/bin/env python

library(stringr)
library(zoo)
library(readr)
library(dplyr)
library(cowplot)
library(purrr)

read_gc_table <- function(filename, datapath){

    samplename = str_split(filename,'\\.')[[1]][1]
    df <- datapath %>%
		str_c(filename, sep='/') %>%
		read_tsv(skip = 6) %>%
        mutate(samplename = samplename) %>%
        mutate(normalize_windows = WINDOWS/sum(WINDOWS, na.rm=T)) %>%
        mutate(subsampled = ifelse(grepl('subsample',filename),'Yes','No'))
    return(df)
}


project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
figure_path  <- str_c(project_path, '/figures')
figurename <- str_c(figure_path, '/gc_plot.pdf')
table_names <- list.files(path = picard_path, pattern = 'gc_metrics')
df <- table_names %>%
	map(read_gc_table, picard_path) %>%
	reduce(rbind) %>%
    filter(samplename != 'SRR733099') %>%
#    filter(!grepl('clustered|sim',samplename)) %>%
#    filter(grepl('kq|q5',samplename)) %>%
    filter(!grepl('X',samplename) ) %>%
    tbl_df

windows_df <- df %>% 
	filter(samplename == unique(.$samplename)[1]) %>%
	mutate(rol_window = normalize_windows * 10) %>%
	mutate(roll_mean_window = rollmean(rol_window, k = 10,fill=0, align='center'))

gc_p <- ggplot(data = df, aes(x = GC, y = NORMALIZED_COVERAGE)) +
	geom_line(size = 1.3, aes(color = samplename)) +
	geom_hline(yintercept = 1, linetype = 2, alpha = 0.9) +
	geom_bar(data = windows_df, aes(x = GC, y = rol_window), 
			 stat='identity', fill='salmon', alpha = 1)  +
	theme(text = element_text(size = 20)) +
	theme(axis.text = element_text(size = 18)) +
#    theme(legend.position =  'None') +
	labs(x = 'GC %', y = 'Normalized Coverage', color = ' ')+
    ylim(0,4)
#ggsave(gc_p, file= figurename)
#message( 'Saved: ' ,figurename)
