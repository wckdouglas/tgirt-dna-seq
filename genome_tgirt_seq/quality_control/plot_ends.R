#!usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)


make_prep <- function(x){
    ifelse(grepl('nextera',x),'Nextera XT','TGIRT-seq')
}

read_files <- function(filename){
    sample_name = str_replace(filename,'.csv','')
    d <- datapath %>%
        str_c(filename,sep='/') %>%
        read_csv() %>%
        select(positions,read_end,base,base_count,base_fraction) %>%
        mutate(filename = sample_name)
    return(d)
}

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/fragment_ends'
files <- list.files(path = datapath, pattern = '.csv')
df <- files %>%
    map(read_files) %>%
    reduce(rbind) %>%
    filter(!grepl('subsampled|sim|Ecoli|q5|phus|SRR',filename)) %>%
    filter(grepl('nextera|NEB|kq|kh|UMI',filename)) %>%
    filter(grepl('MarkDuplicate',filename)) %>%
    filter(grepl('nextera|clustered|umi2id',filename)) %>%
    mutate(prep = case_when(grepl('nextera',.$filename) ~ 'Nextera XT',
                            grepl('pb',.$filename) ~ 'Pacbio',
                            grepl('sim',.$filename) ~ 'Covaris Sim',
                            grepl('SRR',.$filename) ~ 'Covaris SRR',
                            grepl('UMI',.$filename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('NEB',.$filename) ~ 'TGIRT-seq Fragmentase')) %>%
    mutate(prep = ifelse(is.na(prep),'TGIRT-seq Covaris',prep))  %>%
    mutate(read_end = ifelse(read_end == "5'", 'Read 1', 'Read 2')) %>%
    mutate(read_end = factor(read_end, levels=c("Read 1","Read 2"))) %>%
    mutate(actual_positions = ifelse(read_end == "Read 2", positions-20, positions)) %>%
    mutate(bit = -log2(base_fraction)) %>%
    mutate(read_end = as.character(read_end))

colors <- c('light sky blue','salmon','green','yellow')
p <- ggplot(data = df, aes(x = actual_positions, 
                           color = prep, 
                           group=filename, 
                           y = base_fraction)) +
    geom_line(size = 1.3, alpha=0.6) +
    facet_grid(base~read_end, scale ='free_x') +
    labs(x = 'Position Relative to Read ends',y='Fraction of Reads',color=' ') +
    panel_border() +
    scale_color_manual(values = colors)+
    theme(strip.text.x = element_text(size = 20, face='bold')) +
    theme(strip.text.y = element_text(size = 20, face='bold', angle = 0)) +
    theme(axis.title = element_text(size = 20, face='bold')) +
    theme(axis.text = element_text(size = 18, face='bold'))  +
    theme(legend.position = c(0.65,0.45))+
    theme(legend.text = element_text(size = 18, face='bold'))+
    theme(legend.key.size=unit(8,'mm'))
figurename <- str_c(datapath , '/end_bias_plot.pdf')
ggsave(p , file = figurename, height = 8, width = 8)
message('Plotted: ', figurename)

cleave_p <- ggplot(data = df, 
                   aes(x = actual_positions, 
                        group = interaction(filename,base),
                        y = base_fraction))+
    geom_line(alpha=0.8, size = 1.3, aes(color=base))+  
    facet_grid(prep~.)+
    scale_color_manual(values =  RColorBrewer::brewer.pal(9,'Pastel1'))+
    theme(strip.text.x = element_text(size = 20, face='bold')) +
    theme(axis.title = element_text(size = 20, face='bold')) +
    theme(axis.text = element_text(size = 18, face='bold'))  +
    theme(legend.text = element_text(size = 18, face='bold'))+
    theme(legend.key.size=unit(8,'mm')) +
    labs(x = 'Position relative to fragment ends', 
         y = 'Fraction of Base', color = ' ')
figurename <- str_c(datapath, '/cleavage_pattern.pdf')
ggsave(p, file=figurename, height=10,width=10)
message('Plotted: ', figurename)

bit_df <- df %>%
    group_by(actual_positions, read_end, filename, prep) %>%
    summarize(entropy = sum(base_fraction * bit)) %>%
    ungroup()

en_p <- ggplot(data = bit_df, aes(x = actual_positions, y = entropy, 
                      group = filename, color = prep)) + 
    geom_line() +
    facet_grid(.~read_end, scale= 'free_x') +
    geom_line(size = 1, alpha=0.6) +
    labs(x = 'Relative position to fragment ends',y='Entropy (Information)',color=' ') +
    panel_border() +
    scale_color_manual(values = colors) +
    theme(strip.text.x = element_text(size = 20, face='bold')) +
    theme(strip.text.y = element_text(size = 20, face='bold', angle = 0)) +
    theme(axis.title = element_text(size = 20, face='bold')) +
    theme(axis.text = element_text(size = 18, face='bold'))  +
    theme(legend.position = c(0.65,0.45))+
    theme(legend.text = element_text(size = 18, face='bold'))+
    theme(legend.key.size=unit(8,'mm'))+
    ylim(0,2)
figurename <- str_c(datapath, '/ends_entropy.pdf')
ggsave(p, file=figurename, height=10,width=10)
message('Plotted: ', figurename)