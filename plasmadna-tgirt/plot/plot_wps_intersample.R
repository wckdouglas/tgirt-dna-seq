#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(extrafont)
#library(FBN)
loadfonts()

datapath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures'
df <- datapath %>%
    str_c('predictedNucleosomeDistance.tsv',sep='/') %>%
    read_tsv() %>%
    group_by(distance) %>%
    mutate(count = sum(count))


xlim <- 720
internucleosome <- 175
text_x <- 480
text_y <- 1.8
label_xend <- 182
label_yend <- 0.8
label_x <- 360
label_y <- 1.5

#xlim <- 420
nucleo_p <- ggplot(data = df, aes(x = distance, weights=count/1000000)) +
    geom_histogram(fill='black', color = 'black',binwidth = 6, alpha=0.8) +
    scale_x_continuous(breaks = seq(-xlim,xlim,120),limits=c(-xlim,xlim)) +
    theme(text = element_text(size=30, family='Arial', face='plain')) +
    theme(axis.text.x = element_text(size=30,face='plain', family='Arial',angle=50, hjust=0.5, vjust=0.5)) +
    theme(axis.text.y = element_text(size=30,face='plain', family='Arial')) +
    labs(x = 'Inter-nucleosome distance between\ndifferent male individual analyzed by\nTGIRT-seq and ssDNA-seq\n(bp)', 
         y = expression(Peak~count~(x10^{6}))) +
    annotate('segment',xend = label_xend, yend = label_yend, 
             x = label_x, y = label_y, arrow=grid::arrow(), size = 1) +
    annotate('text', x = text_x, y = text_y, label = '180 bp', size = 10) +
    annotate('segment',xend = -label_xend, yend = label_yend, 
             x = -label_x, y = label_y, arrow=grid::arrow(), size = 1) +
    annotate('text', x = -text_x ,y = text_y, label = '-180 bp', size = 10) 

    
label <- expression(paste('x10'^{4}))
nucleo_p <- ggdraw(nucleo_p) 
figurename <- str_c(datapath, '/predictedNucleosomeDistance.pdf')
ggsave(nucleo_p, file=figurename, width = 11,height = 10)
message('Plotted: ', figurename)


distance_probability <-df %>% 
    ungroup() %>%
    mutate(filtered = case_when(
                        abs(.$distance) <= 50 ~ 'one nucleosome',
                        abs(.$distance) > 50 & abs(.$distance <= 190) ~ 'two nucleosome',
                        TRUE ~ 'far away'))  %>% 
    group_by(filtered) %>%
    summarize(total = sum(count)) %>%
    ungroup() %>%
    mutate(fraction = total/sum(total))
