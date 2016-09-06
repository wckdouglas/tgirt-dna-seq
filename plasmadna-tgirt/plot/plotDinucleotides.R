#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)

rename <- function(x){
    y = ifelse(grepl('PD',x),'TGIRT-ssDNA-seq','ssDNA-seq')
    return(y)
}

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA'
data_path <- str_c(project_path,'/nucleotidesAnaylsis/dinucleotides')
figurepath <- str_c(project_path,'/figures')
df <- data_path %>%
    str_c('/dinucleotides.tsv') %>%
    read_tsv()  %>%
    filter(grepl('52|PD',samplename)) %>%
    mutate( name = rename(samplename)) %>%
    filter(position > -120) %>%
    filter(position < 120) %>%
    filter(lenType == 'long (165-170bp)') %>%
    mutate(dinucleotide_type = str_replace_all(dinucleotide_type,"\\|",'/'))
dinucleotide_p <- ggplot(data = df, aes(x = position, y = deviation, color = dinucleotide_type)) +
	    geom_line(size=1.2) +
	    facet_grid(.~name) +
	    theme(axis.text.x = element_text(size=35,angle = 45, hjust= 1,face='plain',family = 'Arial')) +
	    theme(axis.text.y = element_text(size=35, face='plain',family = 'Arial')) +
	    theme(text = element_text(size = 35, face='bold',family = 'Arial')) +
        theme(legend.text = element_text(size = 25, face='plain',family = 'Arial')) + 
        theme(legend.key = element_rect(size = 25)) + 
	    scale_x_continuous(breaks= seq(min(df$position)-1,max(df$position),20)) +
        scale_color_manual(values = c('gold4','blue')) +
	    labs(color = ' ', x = 'Positions relative to center of fragments', y ='Normalized Count') +
        theme(legend.position = c(0.45,0.9))+
        theme(legend.key.height=unit(1.5,"line"))
figurename <- str_c(figurepath, '/dinucleotide.pdf')
ggsave(dinucleotide_p, file= figurename,width=20,height = 8)
message('plotted ', figurename)

	


