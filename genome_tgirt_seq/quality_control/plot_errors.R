#!/usr/bin/env Rscript

library(readr)
library(purrr)
library(cowplot)
library(dplyr)
library(stringr)
library(tidyr)
library(extrafont)
loadfonts()

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
    ifelse(grepl('nextera',x),'Nextera-XT',
           ifelse(grepl('clustered',x),'TGIRT-seq clustered',
                  ifelse(grepl('^K12',x),'TGIRT-seq','TruSeq')))
}

make_prep <- function(x){
    ifelse(grepl('nextera',x),'Nextera-XT','TGIRT-seq')
}

read_data <- function(data_path, annotation){
    tables <- list.files(path = data_path, pattern = '.alignment.metrics', full.names = T)
    tablenames <- data_frame( samplename= tables) %>%
        filter(grepl('nextera|umi2id|cluster',samplename)) %>%
        filter(grepl('^75|UMI|kh|kq|NEB', basename(samplename))) %>%
        tbl_df()
    df <- tablenames$samplename %>%
        map_df(read_table_file) %>%
        mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera-XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('SRR',.$samplename) ~ 'Covaris SRR',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
        mutate(clustered = ifelse(grepl('clustered',samplename),' error-corrected','')) %>%
        mutate(clustered= ifelse(grepl('family',samplename),'error-corrected > 1 member',clustered)) %>%
        mutate(prep = str_c(prep, clustered)) %>%
        select(-clustered) %>%
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
label_list <- c('WGS','Homopolymer length < 4')
df <- map2_df(path_list, 
              label_list, 
              read_data) %>%
    filter(read == 'All Reads')  %>%
    mutate(data_label = factor(data_label, levels = label_list))

base_error <- df %>% 
    filter(prep == 'Nextera-XT') %>%
    group_by(prep, data_label, var) %>%
    summarize(base_error = mean(fraction)) %>%
    ungroup() %>%
    select(-prep) %>%
    tbl_df

plot_df <-  df %>% 
    mutate(prep = case_when(
                            grepl('member',.$prep) ~ 'TGRIT-seq\n(Error-corrected)',
                            grepl('error',.$prep) ~ 'TGRIT-seq (Error-corrected all)',
                            grepl('13N',.$prep) ~ 'TGIRT-seq',
                            grepl('Nextera',.$prep) ~ 'Nextera-XT')) %>%
    inner_join(base_error)  %>%
    mutate(fold = fraction/base_error) %>% 
    mutate(var = ifelse(var == 'Indel Rate','Indel rate','Substitution\nrate')) %>%
    mutate(var = factor(var, levels = rev(unique(var)))) %>%
    filter(prep != 'TGRIT-seq (Error-corrected all)') %>%
#    filter(grepl('Error|Nexter', prep)) %>%
    tbl_df

colors <- c('salmon','black','green4','orange')
error_p <- ggplot(data = plot_df,
            aes(x = prep, y = fraction)) +
    #geom_bar(stat='identity', aes(fill = prep, color= prep), alpha=0.7)+
    geom_jitter(aes(fill = prep, color= prep), alpha=0.7, size = 4)+
    #geom_errorbar(aes(ymin = average_error - sd_error, 
    #                  ymax = average_error + sd_error),
    #              width = 0.25) +
    labs(x = '', y = 'Error rate', fill= ' ', color = ' ') +
    theme(text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(strip.text = element_text(size = 20, family='Arial'))+
    theme(axis.text.x = element_blank()) +
    scale_color_manual(values = colors, guide = guide_legend(ncol=1, keywidth = 2.5, keyheight=1.7)) +
    #scale_fill_manual(values = c('light sky blue','salmon','green4','orange'), guide = guide_legend(ncol=3)) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(legend.text = element_text(size = 20)) +
    #theme(legend.key.size = unit(9, 'mm')) +
    #theme(legend.position = c(0.5, -0.03)) + #poster
    theme(legend.position = c(0.7,0.4)) + #paper
    panel_border() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
    facet_grid(var~data_label, scale = 'free_y')  
label1 <- expression(paste('x10'^{'-3'}))
label2 <- expression(paste('x10'^{'-5'}))
source('~/R/legend_to_color.R')
error_fig <- ggdraw(coloring_legend_text_match(error_p, colors)) +
    theme(plot.margin = grid::unit(c(0,0,5,0),'mm')) +
#    draw_label(label1, x = 0.25, y = 0.897, size = 20, fontface ='bold') + #poster
#    draw_label(label2, x = 0.25, y = 0.45, size = 20, fontface ='bold')  +#poster
    draw_label(label1, x = 0.17, y = 0.87, size = 20, fontface ='bold') + #paper
    draw_label(label2, x = 0.17, y = 0.45, size = 20, fontface ='bold')  #paper
figure_name <- str_c(figure_path,'/genome_errors.pdf')
ggsave(error_fig, file = figure_name, height = 9, width = 13) #poster
message('Plotted: ',figure_name)
