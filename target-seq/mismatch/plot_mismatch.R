#!/usr/bin/env Rscript

library(cowplot)
library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)

table_dir <- '/stor/work/Lambowitz/cdw2854/target-seq/base_tables'
tables <- list.files(path = table_dir, pattern = '.tsv')

read_file <- function(tablename){
    samplename <- str_replace(tablename, '.tsv','')
    df <- table_dir %>%
        str_c(tablename, sep='/') %>%
        read_tsv() %>%
        mutate(samplename = samplename) %>%
        filter(strand == '+') %>%
        filter(start > 174) %>%
        filter(start < 370) %>%
        filter(chrom == 'ENST00000621411')  %>%
#        filter(chrom %in% c("ENST00000613854","ENST00000612966","ENST00000621411",
#                            "ENST00000616365","ENST00000359303","ENST00000356476",
#                            "ENST00000614911","ENST00000618052")) %>%
        tbl_df
    return (df)
}

transition <- c('A-to-G','G-to-A','C-to-T','T-to-C')

assignTemplate <- function(x){
    ifelse(grepl('^H',x),'RNA','DNA')
}

assignFilter <- function(x){
    ifelse(grepl('filtered',x),'>= 3 cluster','No filter')
}

figurename <- str_c(table_dir, '/mismatch.pdf')
merge_df <- tables %>%
    map(read_file) %>%
    reduce(rbind)  %>%
    filter(grepl('bayesian',samplename)) %>%
    filter(grepl('filter',samplename)) %>%
    tbl_df

clean_mismatch_df <- merge_df %>%
    select(ref_base, A,C,T,G, samplename, start) %>%
    gather(base, count, -samplename, -ref_base, -start)  %>%
    tbl_df

mismatch_df <- clean_mismatch_df %>%
    group_by(samplename, base, ref_base) %>%
    summarize(count = sum(count)) %>%
    ungroup() %>% 
    mutate(mutation = str_c(ref_base, '-to-', base)) %>%
    group_by(samplename) %>%
    do(data_frame(
        count = .$count/sum(.$count),
        mutations = .$mutation,
        ref = .$ref_base,
        read = .$base
    )) %>%
    ungroup() %>%
    filter(ref!=read) %>%
    mutate(filter_data = assignFilter(samplename)) %>%
    mutate(template = assignTemplate(samplename) ) %>%
#    mutate(samplename = str_c(template, ' (',filter_data, ')')) %>%
    mutate(label_mutation = ifelse(mutations %in% transition, 'Transition','Transversion')) %>%
    tbl_df

    



### coverage
cov_table <-  merge_df %>% 
    select(samplename, coverage,start)  %>%
    mutate(filter_data = assignFilter(samplename)) %>%
    mutate(template = assignTemplate(samplename) ) %>%
    mutate(samplename = str_c(template, ' (',filter_data, ')')) %>%
    mutate(gene_name = 'HIST1H3B') %>%
    tbl_df
    
cov_p <- ggplot(data = cov_table, aes(x=start, y = coverage)) +
    geom_bar(stat='identity',fill='grey') +
    facet_grid(template~gene_name, scale = 'free_y') +
    theme(text = element_text(size = 20, face='bold')) +
    labs(x = ' ', y = 'High Quality Base Coverage')


pos_mismatch <- clean_mismatch_df %>%
    filter(ref_base != base) %>%
    group_by(samplename, start, ref_base) %>%
    summarize(count = sum(count))   %>%
    ungroup() %>%
    mutate(gene_name = 'HIST1H3B') %>%
    mutate(filter_data = assignFilter(samplename)) %>%
    mutate(template = assignTemplate(samplename) ) %>%
    mutate(samplename = str_c(template, ' (',filter_data, ')')) %>%
    tbl_df

pos_p <- ggplot(data = pos_mismatch, aes(x = start, y = count, fill=ref_base)) +
    geom_bar(stat='identity') +
    facet_grid(template~., scale = 'free_y') +
    theme(text = element_text(size = 20, face='bold')) +
    labs(x = 'Position on GOI', y = 'Mismatch count', fill= ' ') +
    theme(legend.position ='top')


error_DF <- merge_df %>%  
    select(-end,-strand, -insertions, - deletions) %>% 
    gather(base, count, -samplename, -chrom, -start, -ref_base, -coverage) %>% 
    filter(ref_base != base) %>% 
    group_by(samplename, start) %>%
    summarize(count = sum(count),
              coverage = unique(coverage)) %>%
    ungroup() %>%
    group_by(samplename) %>%
    summarize(count = sum(count),
              coverage = sum(coverage)) %>%
    ungroup() %>%
    mutate(error_rate = count / coverage) %>%
    tbl_df

mismatch_df <- mismatch_df %>% 
    inner_join(error_DF %>% select(-count)) %>%
#    mutate(template = str_c(template, '\n(','Error rate = ',signif(error_rate,3),')')) %>%
    tbl_df
    

mismatch_p <- ggplot(data = mismatch_df, aes(x = mutations, y = count, fill = label_mutation)) +
    geom_bar(stat='identity') + 
    facet_grid(template~ref, scale='free') +
    labs(x = ' ', y = 'Fractions', fill = ' ') +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust =0.5)) +
    theme(text = element_text(size = 20, face='bold')) +
    panel_border() +
    theme(legend.position = 'bottom')


p <- ggdraw() +
    draw_plot(cov_p, 0, 0.666, 1, 0.333) +
    draw_plot(pos_p, 0.02, 0.333, 0.98, 0.333) +
    draw_plot(mismatch_p, 0, 0, 1, 0.333) 
#    draw_plot_label(label = c('(B)','(C)','(D)'), 
#                    x = c(0,0,0),
#                    y = c(1,0.666,0.333), 
#                    size = 20)
ggsave(p , file = figurename, width = 10, height = 15)
message('Plotted: ',figurename)

    