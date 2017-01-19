#!/usr/bin/evn Rscript

library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(cowplot)
library(parallel)

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/tss_periodicity'
gene_expression_table <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/rna.csv'
cell_lines_table <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/labels.txt'
cell_lines <- read_tsv(cell_lines_table) %>%
    set_names(c('tissue_type','cells','cell_type','description','rname'))  %>%
    select(tissue_type, cells)  %>% 
    tbl_df
ge <- read_csv(gene_expression_table) %>%
    set_names(c('id','name','cells','TPM','unit')) %>%
    select(-unit)  %>%
    group_by(id,name) %>% 
    do(data_frame(
        cells = .$cells, 
        zeros = sum(.$TPM==0), 
        TPM= .$TPM
    )) %>%
    filter(zeros >= 3) %>%
    filter(TPM>0) %>%
    mutate(TPM = log2(TPM)) %>%
    inner_join(cell_lines)


make_cor_df <- function(filename, datapath){
    periodogram <- datapath %>%
        str_c(filename , sep='/') %>%
        read_tsv() %>%
        filter(periodicity < 199, periodicity > 193) %>%
        group_by(id, name, type) %>%
        summarize(intensity = mean(sqrt(intensity))) %>%
        inner_join(ge) %>%
        group_by(cells, tissue_type) %>%
        summarize(correlation = cor(TPM, intensity, method='pearson')) %>%
        ungroup() %>%
        mutate(samplename = str_replace(filename,'.bed','')) %>%
        tbl_df
    return(periodogram)
}


files <- list.files(path = datapath, pattern='.bed', full.names = F)
df <- files %>%
    mclapply(.,make_cor_df, datapath, mc.cores=12) %>%
    reduce(rbind) %>%
    mutate(abs_cor = abs(correlation)) %>%
    group_by(samplename) %>%
    mutate(rank = rank(abs_cor)) %>%
    tbl_df

tissues <- c("Abdominal" ,'Brain',"Breast/Female Reproductive","Lung","Lymphoid", 
            "Myeloid","Sarcoma", "Skin", "Urinary/Male Reproductive",
            "Other","Primary Tissue")
plot_df <- df %>% 
    #filter(grepl('SRR|merge|RNa',samplename)) %>%
    #filter(grepl('rmdup|1022|1203',samplename)) %>%
    #filter(!grepl('_1_S3',samplename)) %>%
    mutate(tissue_type = factor(tissue_type, levels = tissues)) %>%
    mutate(sample_type = ifelse(grepl('005[12]|^P1',samplename),'Healthy','Breast Cancer')) %>%
    mutate(prep_type = ifelse(grepl('SRR',samplename), 'ssDNA-seq','TGIRT-seq')) %>%
    tbl_df

color_palette <- c('red','green','orange','purple','khaki3','khaki2',
                   'brown','pink','slateblue','darkgrey','skyblue')
tissue_plot <- ggplot(data=plot_df, 
                      aes(x=samplename,y=rank,fill=tissue_type)) + 
    geom_tile(width=0.9, height=0.9, color = 'white') +
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))+
    scale_fill_manual(values = color_palette) +
    labs(y=' ',x = ' ', fill = ' ') +
    ylim(max(plot_df$rank)-20,max(plot_df$rank))+
    theme(axis.line = element_blank()) +
    theme(axis.text.y  = element_blank()) +
    theme(axis.ticks = element_blank()) 
t_d <- data_frame(y=c(10,10,1), 
                  x=c(1, 2, 2))
triangle <- ggplot(t_d) +
    geom_polygon(aes(x=x,y=y), fill = 'grey')+
    theme(axis.line = element_blank()) +
    theme(axis.text  = element_blank()) +
    theme(axis.ticks = element_blank()) +
    labs(x= ' ', y = ' ') 
p <- ggdraw()+
    draw_plot(triangle, 0,0.25,0.1,0.75) +
    draw_plot(tissue_plot, 0.05,0,0.95,1)+
    draw_plot_label(str_c('Rank by Correlation (',max(plot_df$rank),' cells/tissues)'),
                    0,0.2,angle=90)
    


pca_rank <- df %>% 
    select(-abs_cor, -correlation) %>% 
    filter(!grepl('sim|PD[34]',samplename))%>% 
    spread(samplename, rank) %>% 
    select(-cells,-tissue_type) %>% 
    prcomp() %>% 
    .$rotation %>% 
    data.frame() %>% 
    tibble::rownames_to_column('samplename') %>% 
    mutate(prep = ifelse(grepl('SRR',samplename),'ssDNA-seq','TGIRT-seq')) %>%
    ggplot(aes(x=PC1, y =PC2, label = samplename)) + 
        geom_text(aes(color = prep))