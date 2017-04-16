#!/usr/bin/evn Rscript

library(readr)
library(dplyr)
library(cowplot)
library(purrr)
source('~/R/marginal_plot.R')

colors = c('red','darkgreen','chocolate','purple4','gold','khaki',
          'brown','pink','steelblue','skyblue','darkgrey')
tissues_order = c("Abdominal" ,'Brain',"Breast/Female Reproductive","Lung","Lymphoid", 
                 "Myeloid","Sarcoma", "Skin", "Urinary/Male Reproductive" ,
                 "Primary Tissue",'Other')
                 

plot_tissue <- function(datapath){
    df <- datapath %>%
        file.path(.,'cor_table.tsv') %>%
        read_tsv() %>%
        mutate(tissue_type = factor(tissue_type, levels = tissues_order)) %>%
        set_names(c('cells','tissue_type','tgirt','ssdna')) %>%
        mutate(slope = tgirt/ssdna)
    
    corr <- signif(cor(df[,3:4], method='spearman')[2,1],2)
    correlation_label <- sprintf("\"Spearman's\" ~ rho == %0.2f", corr)
    
    ylims <- c(0, 0.015)
    xlims <- c(0, 0.04)
    #tissue_cor_p<-ggplot(data=df, aes(x=SRR2130051_rmdup, y = P1022_UMI_merged, color = tissue_type)) +
    #tissue_cor_p<-ggplot(data=df, aes(x=SRR2130051_rmdup, y = P13_mix_UMI_unique, color = tissue_type)) +
    tissue_cor_p<-ggplot(data=df, aes(x=ssdna, y = tgirt, color = tissue_type)) +
        geom_point(size=10)+
        scale_color_manual(values= colors) +
#        ylim(ylims[1], ylims[2]) +
#        xlim(xlims[1],xlims[2]) +
        labs(y=expression(paste("Pearson's"~rho~"(TGIRT-seq)")),
             x=expression(paste("Pearson's"~rho~"(ssDNA-seq)")),
             color = ' ') +
        theme(legend.key.size = unit(2,'line')) +
        #theme(legend.position = c(0.95,0.5)) +
        theme(legend.text = element_text(size=30, face='plain', family = 'Arial')) +
        theme(axis.text = element_text(size=30, face='plain', family = 'Arial')) +
        theme(text = element_text(size=30, face='plain', family = 'Arial')) +
        annotate(geom='text', label = correlation_label, 
                 parse=T, x = 0.01,y=0.013,
                 size =10, family='Arial') +
#        ggrepel::geom_label_repel(data=df %>% top_n(10,tgirt), aes(label = cells),size=10) +
        geom_rug(size=2, alpha=0.5)
    source('~/R/legend_to_color.R')
    tissue_cor_p <- ggdraw(coloring_legend_text_match(tissue_cor_p, colors))
    return(tissue_cor_p)
}

project_path <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/periodicity_'
analytic_type <- expand.grid(c('tss','gene_body'),c('R','scipy')) %>%
    mutate(folder = str_c(project_path,Var1,'/', Var2))  %>%
    mutate(title = str_c(Var1, '_',Var2))
ps <- lapply(analytic_type$folder,plot_tissue)
p <- plot_grid(plotlist=ps, labels = analytic_type$title)
figurepath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures' 
figure_name <- str_c(figurepath, '/fft_plot.pdf')
ggsave(p, height=30, width=30, file = figure_name)
message('plotted', figure_name)
tissue_cor_p <- ps[[3]]

# 
# dens_y <- ggplot(data=df, aes(x = P1022_UMI_merged, fill=tissue_type))+
#     geom_density(alpha=0.5)+
#     scale_fill_manual(values= colors) +
#     xlim(ylims[1], ylims[2])+
#     theme_void() +
#     theme(legend.position='none') +
#     coord_flip()
# dens_x <- ggplot(data=df, aes(x = SRR2130051_rmdup, fill=tissue_type))+
#     geom_density(alpha=0.5)+
#     scale_fill_manual(values= colors) +
#     xlim(xlims[1], xlims[2]) +
#     theme_void() +
#     theme(legend.position='none')
# ggdraw() +
#     draw_plot(tissue_cor_p, 0,0,0.7,0.7) +
#     draw_plot(dens_x, 0,0.7,0.7,0.3)+
#     draw_plot(dens_x, 0.7,0,0.3,0.7)
#     