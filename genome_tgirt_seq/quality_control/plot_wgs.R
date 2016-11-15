#!/usr/bin/env python

library(zoo)
library(stringr)
library(readr)
library(cowplot)
library(purrr)
library(dplyr)

read_wgs_table <- function(filename, datapath){
    samplename <- str_split(filename,'\\.')[[1]][1]
    df <- datapath %>%
		str_c(filename, sep='/') %>%
		read_tsv(skip = 10) %>%
		filter(!is.na(coverage)) %>%
        mutate(rol_count = count/sum(count) * 100) %>%
#		mutate(rol_count = rollmean(count, 10, fill=0, align='center')) %>%
#		mutate(rol_count = rol_count/ sum(rol_count) * 100) %>%
        mutate(samplename = samplename) %>%
        mutate(subsampled = ifelse(grepl('1M',filename),'Yes','No'))
	message('Read ', filename)
	return(df)
}

project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
figure_path <- str_c( project_path, '/figures')
figurename <- str_c( figure_path, '/wgs_plot.pdf')
table_names <- list.files(path = picard_path , pattern = '.wgs.metrics')
table_names <- table_names[grepl('1M',table_names)]
df <- table_names %>%
	map(read_wgs_table, picard_path) %>%
	reduce(rbind)  %>%
    mutate(line_type = 'WGS') %>%
    select(-subsampled) %>%
    filter(!grepl('sim',samplename)) %>%
    filter(!grepl('clustered',samplename)) %>%
    filter(samplename != 'SRR733099')


base_df <- df %>%
	mutate(bases = count * coverage) %>%
    group_by(samplename) %>%
	summarize(bases = sum(bases),
			  count = sum(count)) %>%
	ungroup() %>%
    mutate(theoretical = bases/count)  %>%
    group_by(samplename) %>% 
    do(data_frame(coverage = seq(100000),
                  rol_count = dpois(seq(100000), .$theoretical) * 100))%>%
    ungroup() %>%
    mutate(line_type = 'Theoretical (Poisson)') 
    
plot_df <- rbind(base_df, df %>%select(-baseq_count, -count)) %>%
    mutate(line_type = factor(line_type, levels = c('WGS','Theoretical (Poisson)')))
    
wgs_p <- ggplot() +
	geom_line(data = plot_df, aes(x = coverage, y = rol_count, 
	                              color = samplename, linetype = line_type), 
	          size = 1.3) + 
	xlim(1,200) +
	theme(legend.position = c(0.6, 0.9)) +
    scale_color_discrete(guide = F)+
    scale_linetype_discrete(guide = guide_legend(ncol = 1))+
	theme(text = element_text(size = 20)) +
	theme(axis.text = element_text(size = 18)) +
	labs(x ='Depth of coverage', y = '% of Genome', 
	     color =' ', linetype = ' ') 
#ggsave(wgs_p, file = figurename)
#message('Saved: ', figurename)


