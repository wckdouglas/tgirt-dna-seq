#!/usr/bin/env Rscript

library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(cowplot)
library(tibble)
library(tidyr)
library(broom)

filepath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/tissueContribution'
figurepath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures'
files <- list.files(path=filepath, pattern = '.tsv')


readTable <- function(filename){
    df <- filepath %>%
        str_c(filename,sep='/') %>%
        read_tsv(col_types = 'cnncncnnnnnnnnnnnnnnnnnnn') %>%
        select(coverage) %>%
        setNames(filename)
}

df <- files %>%
    map(readTable) %>%
    reduce(cbind) %>%
    filter(rowSums(.)!=0) %>%
    tbl_df

cor_df <-df %>% 
    setNames(str_replace(names(.),'.tsv','')) %>%
    cor %>% 
    data.frame %>%
    rownames_to_column('sample2') %>% 
    gather(sample1, cor_coef, -sample2) 

p <- ggplot(data = cor_df, aes(x = sample1, y = sample2, fill = cor_coef)) +
    geom_tile() +
    scale_fill_gradient(low='yellow',high='red')

figurename <- str_c(figure_path,'/cor_plot_dhs.pdf')
ggsave(p,file=figurename)

df[df<5 ] = 0
correlation = cor(df$PD_merged.tsv,df$SRR2130052.tsv)
correlation <- signif(correlation,3)

low_color <- 'yellow'
high_color <- 'red'
p <- ggplot(data =df, aes(y=PD_merged.tsv, x=SRR2130052.tsv ))+
    geom_bin2d(bins=40) +
    scale_x_log10(limits=c(5,1000)) +
    scale_y_log10(limits=c(5,1000)) +
    labs(x = 'TGIRT-seq',y = 'ssDNA-seq (ref 1)', titles='TF Binding Sites') +
    annotate('text',x = 10,y=900, label=str_c('r = ',correlation),size=12) +
    geom_abline(intercept=0, slope=1, color='red') +
    theme(text = element_text(size=20, face='bold')) +
    theme(plot.title = element_text(size=25, face='bold'))  +
    scale_fill_gradient(low = low_color, high = high_color) +
    theme(panel.background = element_rect(fill = low_color))
ggsave(p, file= str_c(figure_path,'/scatter_plot_dhs.png'),dpi=700, h = 5000, w = 2000,type = "cairo-png")
