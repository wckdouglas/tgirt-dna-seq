#!/usr/bin/env python

library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(purrr)
library(ineq)
library(forcats)
library(extrafont)
loadfonts()

read_gc_table <- function(filename, datapath){

    samplename = str_replace(filename,'.gc_metrics','')
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
table_names<- table_names[grepl('^75|sim',table_names)]
df <- table_names %>%
	map(read_gc_table, picard_path) %>%
	purrr::reduce(rbind) %>%
    filter(!grepl('Ecoli',samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera-XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
    mutate(prep = ifelse(is.na(prep),'TGIRT-seq Covaris',prep)) %>%
    tbl_df

windows_df <- df %>% 
	filter(samplename == unique(.$samplename)[1]) %>%
	mutate(rol_window = normalize_windows)

gini_df <- df %>% 
    filter(GC>=12, GC<=80) %>%
    group_by(samplename, prep) %>% 
    summarize(gini=ineq(NORMALIZED_COVERAGE,type='Gini')) %>%
    filter(grepl('nextera|_[EF]_|K12_kh|pb|K12_UMI',samplename)) %>%
    filter(!grepl('clustered', samplename)) %>%
    filter(gini < 0.5) %>%
    ungroup() %>%
    tbl_df

tg = gini_df %>% filter(prep == 'TGIRT-seq 13N') %>% .$gini
xt = gini_df %>% filter(prep == 'Nextera-XT') %>% .$gini
tt = t.test(tg, xt)


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
        geom_bar(data = windows_df, aes(x = GC, y = rol_window*10), 
             stat='identity', fill='springgreen1', alpha = 1)  +
        theme(text = element_text(size=30,face='plain',family = 'Arial')) +
        theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
        theme(legend.text = element_text(size=25,face='plain',family = 'Arial')) +
        scale_linetype_manual(guide='none',values = rep(1,8)) +
        labs(x = '% GC', y = 'Normalized coverage', color = ' ')+
        ylim(0,4)+
        scale_y_continuous(sec.axis = sec_axis(trans = ~.*10, name = '% of 100-bp sliding windows'))
    return(p)
}

index_annotation <- gini_df %>% 
    group_by(prep) %>% 
    summarize(
        mean_gini = mean(gini), 
        sd_gini=sd(gini)
    ) %>%
    ungroup() %>%
    mutate(annotation = str_c('Gini: ', signif(mean_gini,2),' ± ',signif(sd_gini,1))) %>%
    select(prep, annotation) %>%
    tbl_df

gc_df <- df %>% 
#    filter(grepl('nextera|^K12_kh|^K12_kq', samplename)) %>%
#    filter(grepl('nextera|clustered',samplename)) %>%
    filter(grepl('nextera|_[EF]_|K12_kh|pb|K12_UMI',samplename)) %>%
    filter(!grepl('^SRR',samplename)) %>%
    filter(!grepl('SRR1536433',samplename)) %>%
    filter(grepl('UMI|nextera',samplename)) %>% #added for poster
    inner_join(index_annotation) %>% # added for poster
 #   mutate(prep = str_c(prep, '(', annotation, ')')) %>% #added for poster
    mutate(prep = str_replace_all(prep,'13N','')) %>% #added for poster
    mutate(prep = factor(prep, levels = rev(unique(prep)))) 
colors <- c('black','salmon')
gc_p <- plot_gc(gc_df) +
        scale_color_manual(values = colors)+
#        theme(legend.position =  c(0.3,0.7))+ #poster
        theme(legend.position =  c(0.4,0.6))+ #paper
        theme(legend.key.height = unit(2,'line'))
figurename <- str_c(figure_path, '/gc_plot.pdf')
source('~/R/legend_to_color.R')
gc_p<-ggdraw(coloring_legend_text_match(gc_p,colors)) +
  annotate('text', x=0.23, y = 0.41, 
           label = 'No bias', size = 7)
ggsave(gc_p, file = figurename , height = 7, width = 9)
message('Plotted: ', figurename)


linearity <- gc_df %>% 
    filter(GC <= 75, GC>=12) %>% 
    group_by(samplename,prep) %>% 
    nest() %>% 
    mutate(model = map(data, ~lm(NORMALIZED_COVERAGE~GC,data=.))) %>% 
    mutate(summary = map(model,broom::glance)) %>% 
    unnest(summary) %>% 
    group_by(prep) %>%
    summarize(r_mean = mean(r.squared),
              r_sd = sd(r.squared)) %>%
    tbl_df


supplemental_df <- df %>%
    filter(grepl('K12_UMI_3|no_bias|13N',samplename)) %>%
    filter(grepl('[0-9]$', samplename)) %>%
    filter(grepl('13N_K12_sim_template_switch|75bp_K12_UMI_3|13N_K12_sim|13N_K12_sim_ligation_only|no_bias',samplename)) %>%
#    filter(grepl('13N_kmer_K12_sim_template_switch|75bp_K12_UMI_1|13N_kmer_K12_sim|13N_K12_sim_ligation_only|no_bias',samplename)) %>%
    mutate(prep = case_when(grepl('no_bias',.$samplename) ~ 'Simulation: No bias',
                            grepl('sim.[0-9]$',.$samplename) ~'Simulation: Reads 1 and 2 bias',
                            grepl('^fragmentase',.$samplename) ~ 'Simulation: Template switching/Fragmentase bias only',
                            grepl('sim_template_switch',.$samplename) ~ 'Simulation: Template switching',#/Covaris bias only',
                            grepl('ligation',.$samplename) ~ 'Simulation: Ligation bias only',
                            TRUE ~ 'Experimental')) %>%
#    filter(!grepl('13N_K12_sim_template_switch|13N_K12_sim.[0-9]',samplename)) %>%
    mutate(prep = str_replace(prep,'Simulation: ',''))

supplement_df <- supplemental_df %>% 
    filter(prep == 'Experimental') %>% 
    select(GC,NORMALIZED_COVERAGE) %>% 
    group_by(GC) %>%
    summarize(experiment = median(NORMALIZED_COVERAGE)) %>%
    inner_join(supplemental_df) %>%
#    filter(grepl('UMI|.2$',samplename)) %>%
    tbl_df

rmse_df <- supplement_df %>% 
    #filter(GC <= 80, GC>=12) %>% 
    group_by(prep) %>% 
    summarize(rmse = sqrt(mean((experiment-NORMALIZED_COVERAGE)^2))) %>%
    ungroup() %>%
    mutate(rmse = signif(rmse, 3)) %>%
    inner_join(supplement_df) %>%
    mutate(prep = str_c(prep, '\n(RMSE: ',rmse,')')) %>%
    mutate(prep = ifelse(grepl('Experimental',prep),'Experimental', prep)) %>%
    mutate(prep = fct_reorder(prep, rmse)) %>%
    tbl_df

colors <- c('black','red','goldenrod4','springgreen4','navyblue','grey72')
supplemental_p <- plot_gc(rmse_df) + 
        scale_color_manual(values = colors) +
        theme(legend.position = c(0.4,0.75)) +
        theme(legend.key.height = unit(4,'line')) +
        theme(legend.text = element_text(size = 23, hjust=0.5)) 
source('~/R/legend_to_color.R')
supplemental_p <-ggdraw(coloring_legend_text_match(supplemental_p, colors))
figurename <- str_c(figure_path, '/supplemental_gc_plot.pdf')
ggsave(supplemental_p, file = figurename , height = 9, width = 13)
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
    mutate(prep = case_when(grepl('Nextera',.$prep) ~ 'Nextera-XT',
                            grepl('13N',.$prep) ~ 'UMI direct ligation',
                            grepl('Cov',.$prep) ~ 'UMI + CATCG')) %>%
    filter(!grepl('CATCG',prep)) %>%
    mutate(prep = ifelse(grepl('UMI',prep),'TGIRT-seq', prep)) %>%
    tbl_df
lonrenz_df_annotation <- lonrenz_df %>%
    group_by(samplename, prep)  %>%
    summarize(gini=ineq(NORMALIZED_COVERAGE,type='Gini'))  %>%
    ungroup() %>%
    group_by(prep) %>% 
    summarize(
        mean_gini = mean(gini), 
        sd_gini=sd(gini)
    ) %>%
    ungroup() %>%
    mutate(annotation = str_c('Gini: ', signif(mean_gini,2),' ± ',signif(sd_gini,1))) %>%
    select(prep, annotation) %>%
    tbl_df 
lonrenz_df <- lonrenz_df %>% 
    group_by(samplename, prep) %>%
    nest() %>% 
    mutate(lc = map(data,~Lc(.$NORMALIZED_COVERAGE))) %>% 
    mutate(l = map(lc, function(x) x$L)) %>% 
    mutate(p = map(lc, function(x) x$p)) %>% 
    unnest(l,p) %>%
    inner_join(lonrenz_df_annotation) %>%
    mutate(prep = str_c(prep, '(', annotation, ')')) %>% #added for poster
    mutate(prep = factor(prep, levels = rev(unique(prep))))  %>%
    tbl_df
    
colors <- c('black','salmon')
lonrenz_curve <- ggplot(data = lonrenz_df, aes(y= l,x = p, group=samplename, color=prep)) + 
    geom_line(alpha=0.5, size = 1.3) + 
    geom_abline(intercept = 0, slope = 1, linetype=2)+
    scale_x_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25) * (80-12) + 12)+
    theme(text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(axis.text = element_text(size=30,face='plain',family = 'Arial')) +
    theme(legend.text = element_text(size = 25))  +
    theme(legend.key.height =unit(2,'line')) +
    theme(legend.position = c(0.4,0.8))+
    scale_color_manual(values = colors, guide = guide_legend(ncol = 1))+
    labs(x = '% GC', y = 'Cumulative coverage', color = ' ')#, title = 'Lonrenze Curve')
lonrenz_curve <- ggdraw(coloring_legend_text_match(lonrenz_curve, colors))
figurename <- str_c(figure_path, '/gc_lonrez_curve.pdf')
ggsave(lonrenz_curve, file = figurename , height = 7, width = 7)
message('Plotted: ', figurename)
