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
        mutate(group = cut(periodicity, breaks = seq(120,280,10))) %>%
        group_by(group, id, type, name) %>%
        summarize(intensity = sum(intensity)) %>%
        ungroup() %>%
        inner_join(ge) %>%
        group_by(cells, group) %>%
        summarize(correlation = cor(TPM, intensity, method='pearson')) %>%
        ungroup() %>%
        inner_join(cell_lines) %>%
        mutate(cell_kind = ifelse(cells=='NB-4', 'NB4','Others')) %>%
        mutate(periodicity = as.numeric(str_sub(group,2,4))) %>%
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


p<-ggplot(data = df %>%
              filter(samplename %in% c('PD_1203_merged','SRR2130051')), 
          aes(x=periodicity, y = correlation,
              group = cells,color = cell_kind)) +
    geom_line(alpha=0.5) +
   # theme(legend.position = 'none')+
    scale_color_manual(values = c('black','grey')) +
    scale_x_continuous(breaks = seq(160,220,10),limits = c(160,220)) +
    facet_wrap(~samplename)

