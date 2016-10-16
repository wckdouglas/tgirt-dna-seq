#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(cowplot)
library(dplyr)

rename_sample <- function(x){
    ifelse(grepl('^D_',x),'HIST1H3B Plasmid',
           ifelse(grepl('^JR',x),'Jurkat RNA',
                  ifelse(grepl('^JD',x),'Jurkat Genome','Heart RNA')))
}

datapath <- '/stor/work/Lambowitz/cdw2854/target-seq/removed_clipped'
df <- datapath %>%
    str_c('count.tsv',sep='/') %>%
    read_tsv() %>%
    group_by(histone,samplename) %>%
    summarize(normalizedCount = sum(normalizedCount) * 100) %>%
    mutate(template = rename_sample(samplename)) %>%
    mutate(template_type = ifelse(grepl('RNA',template),'Total RNA','DNA'))

p <- ggplot(data = df, aes( x = template, y = normalizedCount, fill=histone, color=histone)) +
    geom_bar(stat='identity') +
    facet_grid(.~template_type, scale='free_x')+
    labs( x = ' ', y = '% of Reads', fill= ' ', color = ' ') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20)) +
    theme(text = element_text(size = 20, face='bold')) +
    theme(axis.text.y = element_text(size = 20, face='bold')) 
figurename <- str_c(datapath, 'gene_count.pdf',sep='/')
ggsave(p, file= figurename ,width = 9, height=7)
message('plotted: ', figurename)
    
