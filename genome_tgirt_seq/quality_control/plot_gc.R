#!/usr/bin/env python

library(stringr)
library(readr)
library(dplyr)
library(cowplot)
library(purrr)
library(ineq)

read_gc_table <- function(filename, datapath){

    samplename = str_split(filename,'\\.')[[1]][1]
    df <- datapath %>%
		str_c(filename, sep='/') %>%
		read_tsv(skip = 6) %>%
        mutate(samplename = samplename) %>%
        mutate(normalize_windows = WINDOWS/sum(WINDOWS, na.rm=T)) %>%
        mutate(subsampled = ifelse(grepl('subsample',filename),'Yes','No'))
    return(df)
}



project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
figure_path  <- str_c(project_path, '/figures')
table_names <- list.files(path = picard_path, pattern = 'gc_metrics')
table_names<- table_names[grepl('^75',table_names)]
df <- table_names %>%
	map(read_gc_table, picard_path) %>%
	reduce(rbind) %>%
    filter(!grepl('Ecoli',samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
    mutate(prep = ifelse(is.na(prep),'TGIRT-seq Covaris',prep)) %>%
    tbl_df

windows_df <- df %>% 
	filter(samplename == unique(.$samplename)[1]) %>%
	mutate(rol_window = normalize_windows * 10) %>%
	mutate(roll_mean_window = zoo::rollmean(rol_window, k = 10,fill=0, align='center'))

gini_df <- df %>% 
    filter(GC>=12, GC<=80) %>%
    group_by(samplename, prep) %>% 
    summarize(gini=ineq(NORMALIZED_COVERAGE,type='Gini')) %>%
    filter(grepl('nextera|_[EF]_|K12_kh|pb|K12_UMI',samplename)) %>%
    filter(!grepl('clustered', samplename)) %>%
    filter(gini < 0.5) %>%
    ungroup() %>%
    tbl_df

gini_p <- ggplot(data=gini_df, aes(x = prep, y=gini, color=prep)) +
    geom_point(size = 4) +
    #geom_violin(size = 2) +
    labs(x = ' ', y = 'Gini Coefficient', color = ' ') +
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    theme(text = element_text(size = 25, face='bold')) +
    theme(axis.text = element_text(size = 25, face='bold')) +
    theme(legend.key.height = unit(2,'line')) +
    theme(legend.position='bottom')+
    guides(col = guide_legend(ncol=2))
figurename <- str_c(figure_path, '/gc_gini_plot.pdf')
ggsave(gini_p, file = figurename , height = 7, width = 7)
message('Plotted: ', figurename)



plot_gc <-function(df){
    p <- ggplot(data = df, aes(x = GC, y = NORMALIZED_COVERAGE)) +
        geom_line(size = 1.3, aes(color = prep, group = samplename), alpha=0.7) +
        geom_hline(yintercept = 1, linetype = 2, alpha = 0.9) +
        geom_bar(data = windows_df, aes(x = GC, y = rol_window), 
             stat='identity', fill='springgreen1', alpha = 1)  +
        theme(text = element_text(size = 25, face='bold')) +
        theme(axis.text = element_text(size = 25, face='bold')) +
        scale_linetype_manual(guide='none',values = rep(1,8)) +
        labs(x = 'GC %', y = 'Normalized coverage', color = ' ')+
        ylim(0,4)
    return(p)
}

index_annotation <- gini_df %>% 
    group_by(prep) %>% 
    summarize(
        mean_gini = mean(gini), 
        sd_gini=sd(gini)
    ) %>%
    ungroup() %>%
    mutate(annotation = str_c('Gini: ', signif(mean_gini,3),'±',signif(sd_gini,3))) %>%
    select(prep, annotation) %>%
    tbl_df

gc_p <- df %>% 
#    filter(grepl('nextera|^K12_kh|^K12_kq', samplename)) %>%
#    filter(grepl('nextera|clustered',samplename)) %>%
    filter(grepl('nextera|_[EF]_|K12_kh|pb|K12_UMI',samplename)) %>%
    filter(!grepl('^SRR',samplename)) %>%
    filter(!grepl('SRR1536433',samplename)) %>%
    filter(grepl('UMI|nextera',samplename)) %>% #added for poster
    inner_join(index_annotation) %>% # added for poster
    mutate(prep = str_c(prep, '(', annotation, ')')) %>% #added for poster
    mutate(prep = str_replace_all(prep,'13N','')) %>% #added for poster
    plot_gc() +
        scale_color_manual(values = c('light sky blue','salmon'))+
#        theme(legend.position =  c(0.3,0.7))+ #poster
        theme(legend.position =  c(0.55,0.9))+ #paper
        theme(legend.key.height = unit(2,'line'))
figurename <- str_c(figure_path, '/gc_plot.pdf')
source('~/R/legend_to_color.R')
#gc_p<-ggdraw(coloring_legend_text(gc_p))
ggsave(gc_p, file = figurename , height = 7, width = 9)
message('Plotted: ', figurename)

rename_sim <- function(x){
}

supplemental_df <- df %>%
    filter(grepl('K12_kh|sim',samplename)) %>%
    mutate(prep = case_when(grepl('no_bias',.$samplename) ~ 'Simulation: no bias',
                            grepl('sim$',.$samplename) ~'Simulation: Reads 1 and 2 bias',
                            grepl('sim_template_switch',.$samplename) ~ 'Simulation: Read 1 bias only',
                            grepl('ligation',.$samplename) ~ 'Simluation: Read 2 bias only')) %>%
    mutate(prep = ifelse(is.na(prep),'Experimental',prep)) %>%
    mutate(prep = factor(prep, levels = c('Experimental',
                                           'Simulation: no bias',
                                           'Simulation: Read 1 bias only',
                                           'Simluation: Read 2 bias only',
                                           'Simulation: Reads 1 and 2 bias')))
supplemental_p <- plot_gc(supplemental_df) + 
        scale_color_manual(values = c('grey','skyblue1','salmon','bisque4','red')) +
        theme(legend.position = c(0.3,0.8))
figurename <- str_c(figure_path, '/supplemental_gc_plot.pdf')
ggsave(supplemental_p, file = figurename , height = 8, width = 10)
message('Plotted: ', figurename)
sim_gini <- supplemental_df %>%
    filter(GC>10, GC<75) %>%
    group_by(samplename, prep) %>% 
    summarize(gini=ineq(NORMALIZED_COVERAGE,type='Gini')) %>%
    filter(gini < 0.7) %>%
    ungroup() %>%
    tbl_df

csdf <- supplemental_df %>% 
    filter(GC>10, GC<75) %>%
    arrange(GC) %>%
    group_by(samplename, prep) %>% 
    do(data_frame(
        gc= .$GC, 
        cc = cumsum(.$NORMALIZED_COVERAGE)
    )) %>%
    ungroup()
cs_p <- ggplot(csdf, aes(x=gc, y = cc, color=samplename)) + 
    geom_line() + 
    geom_abline(intercept=-10, slope=1)




lonrenz_df <-df %>% 
    filter(GC>=12, GC<=80) %>% 
    filter(grepl('nextera|_[EF]_|K12_kh|pb|K12_UMI',samplename)) %>%
    filter(!grepl('clustered', samplename)) %>%
    filter(!grepl('SRR1536433',samplename)) %>%
    filter(grepl('nextera|UMI|K12_[kh]',samplename)) %>%
    mutate(prep = case_when(grepl('Nextera',.$prep) ~ 'Nextera XT',
                            grepl('13N',.$prep) ~ 'UMI direct ligation',
                            grepl('Cov',.$prep) ~ 'UMI + CATCG')) %>%
    group_by(samplename, prep) %>% 
    nest() %>% 
    mutate(lc = map(data,~Lc(.$NORMALIZED_COVERAGE))) %>% 
    mutate(l = map(lc, function(x) x$L)) %>% 
    mutate(p = map(lc, function(x) x$p)) %>% 
    unnest(l,p)
lonrenz_curve <- ggplot(data = lonrenz_df, aes(y= l,x = p, group=samplename, color=prep)) + 
    geom_line(alpha=0.5) + 
    geom_abline(intercept = 0, slope = 1)+
    scale_x_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25) * (80-12) + 12)+
    labs(x = '% of GC', y = 'Cumulative coverage', color = ' ')#, title = 'Lonrenze Curve')
figurename <- str_c(figure_path, '/gc_lonrez_curve.pdf')
ggsave(lonrenz_curve, file = figurename , height = 7, width = 7)
message('Plotted: ', figurename)
