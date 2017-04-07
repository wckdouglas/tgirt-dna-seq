#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(cowplot)
library(ggforce)

read_file <- function(filename){
    read_csv(filename) %>%
    mutate(samplename = str_split(basename(filename), '\\.')[[1]][1]) %>%
    return
}

datapath <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/fragment_ends/'
fig_name <- str_c(datapath, '/entropy_plot.pdf')
base_fig_name <- str_c(datapath, '/base_plot.pdf')
combined_fig_name <- str_c(datapath, '/zoomed_entropy_plot.pdf')
fig_width <- 10 
fig_height <- 7
df <- list.files(path = datapath, 
                        pattern = 'MarkDuplicate.csv$',
                        full.names = T) %>%
    .[grepl('kh|kq|7N|UMI',.)] %>%
    .[grepl('umi2id', .)] %>%
    map_df(read_file) %>%
    mutate(adjusted_position = ifelse(read_end == "5'", positions+1, positions-20)) %>%
    mutate(entropy = -base_fraction * log2(base_fraction)) %>%
    mutate(method = case_when(
                    grepl('7N', .$samplename) ~ '7N',
                    grepl('UMI', .$samplename) ~ 'UMI direct ligation',
                    grepl('K12_k', .$samplename) ~ "5' CGATG + UMI")) %>%
    mutate(read_end = ifelse(read_end == "5'", 'Read 1', 'Read 2')) %>%
    filter(!grepl('7N',method)) %>%
    tbl_df()

plot_base_dist <- function(meth){
    p <- ggplot(data = df %>% filter(method == meth), 
                aes(x = adjusted_position, y = base_fraction, color = base))+
        geom_line()+
        facet_grid(samplename~read_end, scale='free_x') +
        labs(x = 'Position', 
               y = 'Base fraction', 
               color = ' ', 
               title = meth)+
        theme(strip.text.y = element_blank())+
        theme(text = element_text(face='bold', size = 25))
}
ps <- lapply(df$method%>%unique(), plot_base_dist)
base_p <- plot_grid(plotlist= ps, ncol = 3)
ggsave(base_p, file = base_fig_name, 
       width=14, height = 10)
message('Plotted ', base_fig_name)


df <- df %>% 
    group_by(samplename, read_end, adjusted_position, method) %>%
    summarize(entropy = sum(entropy)) %>%
    ungroup()
    tbl_df
    
    
en_p <- ggplot(data = df, aes(x = adjusted_position, 
                      y = entropy, group=samplename))+
    geom_line(alpha=0.4, size = 1.5, aes(color = method))+
    facet_grid(.~read_end, scale='free_x')+
    labs(x = 'Positions', y = 'Entropy (bits)', color = ' ') +
    ylim(0,2)+
    theme(text = element_text(face='bold', size = 25)) +
    theme(legend.position = c(0.7,0.5))
ggsave(en_p, file = fig_name, 
       width=fig_width, height = fig_height)
message('Plotted ', fig_name)

colors <- c('black','khaki4','green4')
small_en_p <- ggplot(data = df %>% 
                         filter(read_end=='Read 1') %>% 
                         mutate(method = relevel(factor(method),'UMI direct ligation')), 
                     aes(color = method, x = adjusted_position, 
                        y = entropy, group=samplename)) +
    geom_line(size = 1.5, alpha=0.4) +
    facet_zoom(x = adjusted_position < 5)  +
    theme(legend.key.height = unit(2,'line')) +
    geom_vline(aes(xintercept = adjusted_position), linetype=2, alpha=0.3, color = 'grey') +    
    labs(x = 'Positions', y = 'Entropy (bits)\n', color = ' ')  +
    theme(text = element_text(family='Arial', size = 25)) +
    theme(legend.position = c(0.3,0.15)) +
    theme(axis.text = element_text(family ='Arial', size = 25)) +
    scale_color_manual(values = colors) +
    theme(legend.text = element_text(family ='Arial', size = 25)) 
source('~/R/legend_to_color.R')
small_en_p <- ggdraw(coloring_legend_text_match(small_en_p, colors))
ggsave(small_en_p, file = combined_fig_name, 
       width=fig_width, height = fig_height)
message('Plotted ', combined_fig_name)
    