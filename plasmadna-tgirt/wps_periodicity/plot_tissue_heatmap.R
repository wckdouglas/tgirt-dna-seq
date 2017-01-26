#!/usr/bin/evn Rscript

library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(cowplot)
library(parallel)
library(RColorBrewer)


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
        non_zeros = sum(.$TPM!=0), 
        TPM= .$TPM
    )) %>%
    ungroup() %>%
    filter(non_zeros >= 3) %>%
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
    filter(grepl('SRR|UMI',samplename)) %>%
    filter(grepl('rmdup|UMI_merge',samplename)) %>%
    #filter(grepl('rmdup|SRR|merge',samplename)) %>%
    #filter(!grepl('clustered_rmdup|clustered|S[0-9]_rmdup|[345]$',samplename)) %>%
    mutate(tissue_type = factor(tissue_type, levels = tissues)) %>%
    mutate(sample_type = case_when(
            grepl('005[12]|^P1|^PD',.$samplename)~'Healthy',
            grepl('SRR2130004',.$samplename)~'Breast cancer (Invasive/infiltrating ductal)',
            grepl('0011|0032|0045',.$samplename)~'Breast cancer (Invasive/infiltrating lobular)',
            grepl('0043|0033',.$samplename)~'Breast cancer (Ductal carcinoma in situ)'
        )) %>%
    mutate(prep_type = ifelse(grepl('SRR',samplename), 'ssDNA-seq','TGIRT-seq')) %>%
    mutate(name =samplename) %>%
    mutate(samplename = str_c(sample_type,': ', prep_type, '-',samplename)) %>%
    tbl_df

color_palette <- c('red','green','orange','purple','lightgoldenrod2','khaki1',
                   'brown','pink','slateblue','darkgrey','skyblue')
tissue_plot <- ggplot(data=plot_df, 
                      aes(x=samplename,y=rank,fill=tissue_type)) + 
    geom_tile(width=0.9, height=0.9, color = 'white') +
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))+
    scale_fill_manual(values = color_palette) +
    labs(y=' ',x = ' ', fill = ' ') +
    ylim(max(plot_df$rank)-25,max(plot_df$rank))+
    theme(axis.line = element_blank()) +
    theme(axis.text.y  = element_blank()) +
    theme(axis.ticks = element_blank())  +
    theme(legend.text = element_text(size = 25, face='bold'))+
    theme(legend.key.size  = unit(2, 'line'))
t_d <- data_frame(y=c(1,1,10), 
                  x=c(1, 2, 2))
triangle <- ggplot(t_d) +
    geom_polygon(aes(x=x,y=y), fill = 'grey')+
    theme(axis.line = element_blank()) +
    theme(axis.text  = element_blank()) +
    theme(axis.ticks = element_blank()) +
    labs(x= ' ', y = ' ') 
label_coloring <- ggplot(data= plot_df %>% 
                             select(samplename,sample_type,prep_type)  %>%
                             unique() %>%
                             gather(category, fill_type, -samplename))+    
    geom_tile(aes(y = category, x = samplename, fill = fill_type)) +
    theme(axis.line = element_blank()) +
    theme(axis.text  = element_blank()) +
    theme(axis.ticks = element_blank()) +
    labs(y= ' ', x=' ', fill = ' ') +
    theme(legend.position = 'bottom') +
    scale_fill_manual(values= brewer.pal(8,'Dark2')) +
    theme(legend.text = element_text(size = 25, face='bold')) +
    theme(legend.key.size  = unit(2, 'line'))

p <- ggdraw()+
    draw_plot(triangle, 0.02,0.09,0.06,0.7) +
    draw_plot(tissue_plot + theme(axis.text.x = element_blank(),
                                  legend.position = 'top') , 
              0.05,0.07,0.95,0.82)+
    draw_plot(label_coloring, 0.05, 0 ,0.95,0.14) +
    draw_plot_label(str_c('Rank by Correlation (',max(plot_df$rank),' cells/tissues)'),
                    0,0,angle=90, size =30, hjust = -0.4) 
figurename <- str_c(datapath,'tissue_inference.pdf',sep='/')
ggsave(p, file = figurename, width=20,height=15)
message('Plotted: ', figurename)




pca_rank <- df %>% 
    select(-correlation, -abs_cor) %>% 
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

cor_df <- df %>%
    filter(!grepl('rmdup',samplename))%>%
    select(-correlation, - abs_cor) %>% 
    spread(samplename, rank) %>% 
    select(-cells,-tissue_type) %>% 
    cor %>%
    data.frame() 

heatmap3::heatmap3(cor_df, col = viridis::viridis(1000),
         symm=T)
#    tibble::rownames_to_column('sample1') %>%
#    gather(sample2, cor_value,-sample1) %>%
#    ggplot(aes(x=sample1, y = sample2, fill = cor_value)) +
#        geom_tile() + 
#        viridis::scale_fill_viridis()