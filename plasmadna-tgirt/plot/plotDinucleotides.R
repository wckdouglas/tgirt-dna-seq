#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(extrafont)

read_file <- function(tablename){
    
}

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA'
data_path <- str_c(project_path,'/nucleotidesAnaylsis/dinucleotides')
figurepath <- str_c(project_path,'/figures')
filenames <- list.files(path = data_path, 
                        pattern = '.tsv', 
                        full.names = T)

df <- filenames %>%
    #.[grepl('SRR2130051.tsv|P1022_2_S4_umi2id_unique.tsv',.)] %>%
    .[grepl('51_rmdup|P1022_1113_13_1016_mix_unique',.)]%>%
    map_df(read_tsv) %>%
    #filter(grepl('52|P1113_3_S20_umi2id', samplename)) %>%
    mutate( name = ifelse(grepl('SRR',samplename),'ssDNA-seq','TGIRT-seq')) %>%
    filter(position > -120) %>%
    filter(position < 120) %>%
    mutate(dinucleotide_type = str_replace_all(dinucleotide_type,"\\|",'/')) %>%
    mutate(name = factor(name, levels = c('TGIRT-seq','ssDNA-seq')))
dinucleotide_p <- ggplot(data = df, aes(x = position, y = adjusted_signal, color = dinucleotide_type)) +
	    geom_line(size=1.2) +
	    facet_grid(.~name) +
	    theme(axis.text.x = element_text(size=30,angle = 45, hjust= 1,face='plain',family = 'Arial')) +
	    theme(axis.text.y = element_text(size=30, face='plain',family = 'Arial')) +
	    theme(text = element_text(size = 30, family = 'Arial')) +
        theme(legend.text = element_text(size = 30, family = 'Arial')) + 
	    scale_x_continuous(breaks= seq(min(df$position)-1,max(df$position),20)) +
        scale_color_manual(values = c('gold4','blue')) +
	    labs(color = ' ', x = 'Position relative to center of 167-nt fragments (bp)', y ='Normalized count') +
        theme(legend.position = c(0.5,0.9))+
        theme(legend.key.height=unit(2,"line"))
source('~/R/legend_to_color.R')
dinucleotide_p <- ggdraw(coloring_legend_text(dinucleotide_p))
figurename <- str_c(figurepath, '/dinucleotide.pdf')
ggsave(dinucleotide_p, file= figurename,width=20,height = 8)
message('plotted ', figurename)

	


