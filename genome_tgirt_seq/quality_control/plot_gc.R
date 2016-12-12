#!/usr/bin/env python

library(stringr)
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

make_prep <- function(x){
    ifelse(grepl('nextera',x),'Nextera XT',
           ifelse(grepl('pb',x),'Pacbio','TGIRT-seq'))
}


project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
figure_path  <- str_c(project_path, '/figures')
figurename <- str_c(figure_path, '/gc_plot.pdf')
table_names <- list.files(path = picard_path, pattern = 'gc_metrics')
df <- table_names %>%
	map(read_gc_table, picard_path) %>%
	reduce(rbind) %>%
#    dplyr::filter(grepl('nextera|^K12_kh|^K12_kq', samplename)) %>%
    mutate(prep = make_prep(samplename)) %>%
#    filter(!grepl('clustered|sim',samplename)) %>%
#    filter(grepl('kq|q5',samplename)) %>%
#    filter(!grepl('X',samplename) ) %>%
    tbl_df

windows_df <- df %>% 
	filter(samplename == unique(.$samplename)[1]) %>%
	mutate(rol_window = normalize_windows * 10) %>%
	mutate(roll_mean_window = zoo::rollmean(rol_window, k = 10,fill=0, align='center'))


plot_gc <-function(df){
    p <- ggplot(data = df, aes(x = GC, y = NORMALIZED_COVERAGE)) +
        geom_line(size = 1.3, aes(color = prep, group = samplename)) +
        geom_hline(yintercept = 1, linetype = 2, alpha = 0.9) +
        geom_bar(data = windows_df, aes(x = GC, y = rol_window), 
             stat='identity', fill='springgreen1', alpha = 1)  +
        theme(text = element_text(size = 20)) +
        theme(axis.text = element_text(size = 18)) +
        scale_linetype_manual(guide='none',values = rep(1,8)) +
        labs(x = 'GC %', y = 'Normalized Coverage', color = ' ')+
        ylim(0,4)
    return(p)
}

gc_p <- df %>% 
    filter(grepl('nextera|^K12_kh|^K12_kq', samplename)) %>%
    filter(grepl('nextera|clustered',samplename)) %>%
    plot_gc() +
        theme(legend.position =  'None') 

rename_sim <- function(x){
    ifelse(grepl('no_bias',x), 'Simulation: no bias',
           ifelse(grepl('sim$',x),'Simulation: Reads 1 and 2 bias',
                ifelse(grepl('sim_template_switch',x),'Simulation: Read 1 bias only',
                        ifelse(grepl('ligation',x),'Simluation: Read 2 bias only','Experimental'))))
}

supplemental_p <- df %>%
    filter(grepl('K12_kh|sim',samplename)) %>%
    mutate(prep = rename_sim(samplename)) %>%
    mutate(prep = factor(prep, levels = c('Experimental',
                                           'Simulation: no bias',
                                           'Simulation: Read 1 bias only',
                                           'Simluation: Read 2 bias only',
                                           'Simulation: Reads 1 and 2 bias'))) %>%
    plot_gc() + 
        scale_color_manual(values = c('grey','skyblue1','salmon','bisque4','red')) +
        theme(legend.position = c(0.3,0.8))
figurename <- str_c(figure_path, '/supplemental_gc_plot.pdf')
ggsave(supplemental_p, file = figurename , height = 8, width = 10)
message('Plotted: ', figurename)
        