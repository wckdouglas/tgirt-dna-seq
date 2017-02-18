#!/usr/bin/env Rscript

library(readr)
library(purrr)
library(cowplot)
library(dplyr)
library(stringr)
library(tidyr)

rename_enzyme <- function(x){
    if (grepl('kh',x)){
        'KAPA HIFI'
    }else if (grepl('kq',x)) {
        'KAPA QUANT'
    }else if (grepl('ph',x)){
        'Phusion'
    }else if (grepl('q5',x)){
        'Q5'
    }else{
        ' '
    }
}

project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
data_path1 <- str_c( project_path, '/picard_results')
data_path2 <- str_c( project_path, '/picard_results/no_hiINDEL')
figure_path  <- str_c(project_path, '/figures')

read_table_file <- function(tablename){
    samplename <- str_replace(basename(tablename),'.MarkDuplicate.subsampled.alignment_metrics','')
    message('Read ', samplename)
    df <- read_tsv(tablename, skip=6) %>% 
        head(3) %>% 
        gather(var, value, -CATEGORY) %>% 
        spread(CATEGORY,value) %>%
        mutate(samplename = samplename) %>%
        filter(var %in% c('PF_MISMATCH_RATE','PF_INDEL_RATE')) 
    return(df)
}

rename_sample <- function(x){
    ifelse(grepl('nextera',x),'Nextera XT',
           ifelse(grepl('clustered',x),'TGIRT-seq clustered',
                  ifelse(grepl('^K12',x),'TGIRT-seq','TruSeq')))
}

make_prep <- function(x){
    ifelse(grepl('nextera',x),'Nextera XT','TGIRT-seq')
}

read_data <- function(data_path, annotation){
    tablenames <- list.files(path = data_path, pattern = '.alignment.metrics', full.names = T)
    df <- tablenames %>%
        map_df(read_table_file) %>%
        filter(grepl('nextera|umi2id|cluster',samplename)) %>%
        filter(grepl('^75|UMI|kh|kq|NEB',samplename)) %>%
        mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('SRR',.$samplename) ~ 'Covaris SRR',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
        mutate(clustered = ifelse(grepl('clustered',samplename),' error-corrected','')) %>%
        mutate(prep = str_c(prep, clustered)) %>%
        gather(pair, fraction, -samplename, -prep, - var) %>%
        filter(!is.na(fraction)) %>%
        mutate(read = ifelse(grepl('FIRST',pair),'Read 1',ifelse((pair=='PAIR'),'All Reads','Read 2'))) %>%
        mutate(fraction = as.numeric(fraction)) %>%
        #group_by(prep, read, var) %>%
        #summarize(
        #    average_error = mean(fraction),
        #    sd_error = sd(fraction)
        #) %>%
        ungroup() %>%
        mutate(var = str_replace_all(var,'_',' ')) %>%
        mutate(var = str_replace_all(var,'PF ','')) %>%
        mutate(var = str_to_title(var)) %>%
        mutate(fraction = ifelse(grepl('Indel',var), fraction * 1e5, fraction*1e3)) %>%
        #mutate(average_error = ifelse(grepl('Indel',var), average_error * 1e5, average_error*1e3)) %>%
        #mutate(sd_error = ifelse(grepl('Indel',var), sd_error * 1e5, sd_error*1e3)) %>%
        filter(grepl('13N|Nextera', prep)) %>%
        mutate(data_label = annotation) %>%
        tbl_df
}

path_list <- c(data_path1,data_path2)
label_list <- c('WGS','(Homopolymer length < 4)')
df <- map2_df(path_list, 
              label_list, 
              read_data) %>%
    filter(read == 'All Reads')  %>%
    mutate(data_label = factor(data_label, levels = label_list))

base_error <- df %>% 
    filter(prep == 'Nextera XT') %>%
    group_by(prep, data_label) %>%
    summarize(base_error = mean(fraction)) %>%
    ungroup() %>%
    select(-prep)
    tbl_df

plot_df <-  df %>% 
    mutate(prep = case_when(grepl('error',.$prep) ~ 'TGRIT-seq (Error-corrected )',
                            grepl('13N',.$prep) ~ 'TGIRT-seq',
                            grepl('Nextera',.$prep) ~ 'Nextera XT')) %>%
    inner_join(base_error)  %>%
    mutate(fold = fraction/base_error) %>%
    tbl_df

p <- ggplot(data = plot_df,
            aes(x = prep, y = fold)) +
    #geom_bar(stat='identity', aes(fill = prep, color= prep), alpha=0.7)+
    geom_jitter(aes(fill = prep, color= prep), alpha=0.7, size = 5)+
    #geom_errorbar(aes(ymin = average_error - sd_error, 
    #                  ymax = average_error + sd_error),
    #              width = 0.25) +
    labs(x = '', y = 'Error rate relative to Nextera-XT (fold)', fill= ' ', color = ' ') +
    theme(strip.text = element_text(size = 25, face = 'bold'))+
    theme(text = element_text(size = 25, face = 'bold'))+
    theme(axis.text.x = element_blank()) +
    scale_color_manual(values = c('light sky blue','salmon','green4','orange'), guide = guide_legend(ncol=3)) +
    #scale_fill_manual(values = c('light sky blue','salmon','green4','orange'), guide = guide_legend(ncol=3)) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.y = element_text(size = 25, face='bold')) +
    theme(legend.key.size = unit(9, 'mm')) +
    theme(legend.position = c(0.5, -0.01)) + #poster
    #theme(legend.position = 'bottom') + #paper
    panel_border() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
    facet_grid(var~data_label, scale = 'free_y') 
label1 <- expression(paste('x10'^{'-5'}))
label2 <- expression(paste('x10'^{'-3'}))
source('~/R/legend_to_color.R')
#fig <- ggdraw(coloring_legend_text(p)) +
fig <- ggdraw(p) +
#    draw_label(label1, x = 0.25, y = 0.897, size = 20, fontface ='bold') + #poster
#    draw_label(label2, x = 0.25, y = 0.45, size = 20, fontface ='bold')  #poster
    draw_label(label1, x = 0.14, y = 0.897, size = 20, fontface ='bold') + #paper
    draw_label(label2, x = 0.14, y = 0.45, size = 20, fontface ='bold')  #paper
figure_name <- str_c(figure_path,'/genome_errors.pdf')
#ggsave(fig, file = figure_name, height = 7.6, width = 6) #poster
ggsave(fig, file = figure_name, height = 7.6, width = 12)
message('Plotted: ',figure_name)