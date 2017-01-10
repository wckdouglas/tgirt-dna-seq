#!/usr/bin/evn Rscript

library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(cowplot)
library(parallel)

gene_expression <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/rna.csv'
ge <- read_csv(gene_expression) %>%
    set_names(c('id','name','cells','TPM','unit')) %>%
    select(-unit) 
cell_lines <- read_tsv('/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/labels.txt')  %>%
    set_names(c('tissue_type','cells','cell_type','description','rname')) 

make_cor_df <- function(filename, datapath){
    periodogram <- datapath %>%
        str_c(filename , sep='/') %>%
        read_tsv() %>%
        filter(periodicity < 199, periodicity > 193) %>%
        group_by(id, name, type) %>%
        summarize(intensity = mean(intensity)) %>%
        inner_join(ge) %>%
        group_by(cells) %>%
        summarize(correlation = cor(TPM, intensity, method='pearson')) %>%
        ungroup() %>%
        inner_join(cell_lines) %>%
        arrange(abs(correlation)) %>%
        mutate(rank = seq(1:nrow(.))) %>%
        mutate(samplename = str_replace(filename,'.bed','')) %>%
        tbl_df
    return(periodogram)
}

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/tss_periodicity'
files <- list.files(path = datapath, pattern='.bed', full.names = F)
df <- files %>%
    mclapply(.,make_cor_df, datapath, mc.cores=12) %>%
    reduce(rbind) %>%
    tbl_df

tissue_plot <- ggplot(data=tissue_df, 
                      aes(x=samplename,y=rank,fill=tissue_type)) + 
    geom_tile() +
    theme(axis.text.x = element_text(angle=90))+
    ylim(40,60)