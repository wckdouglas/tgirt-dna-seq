#!/usr/bin/env python

library(stringr)
library(readr)
library(cowplot)
library(purrr)
library(dplyr)
library(broom)
library(tidyr)


rename_enzyme <- function(x){
    if (x == 'kh'){
        'KAPA HIFI'
    }else if (x == 'kq') {
        'KAPA QUANT'
    }else if (x == 'ph'){
        'Phusion'
    }else if (x == 'q5'){
        'Q5'
    }else{
        x
    }
}

read_wgs_table <- function(filename, datapath){
    samplename <- str_split(filename,'\\.')[[1]][1]
    df <- datapath %>%
		str_c(filename, sep='/') %>%
		read_tsv(skip = 10) %>%
        filter(coverage<100) %>%
		filter(!is.na(coverage)) %>%
        mutate(density = count/sum(count)) %>%
        mutate(samplename = samplename) %>%
        tbl_df
	message('Read ', filename)
	return(df)
}

project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
figure_path <- str_c( project_path, '/figures')
figurename <- str_c( figure_path, '/wgs_plot.pdf')
table_names <- list.files(path = picard_path , pattern = '.wgs.metrics')
    
df <- table_names %>%
	map(read_wgs_table, picard_path) %>%
    reduce(rbind) %>%
    mutate(line_type = 'WGS') %>%
#    select(-subsampled) %>%
    dplyr::filter(grepl('nextera|^K12_kh|^K12_kq', samplename)) %>%
    dplyr::filter(grepl('nextera|umi2id', samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera XT',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris'))%>%
    tbl_df

base_df <- df %>%
    mutate(bases = count * coverage) %>%
    group_by(samplename,prep) %>%
    summarize(bases = sum(bases),
              count = sum(count)) %>%
    ungroup() %>%
    mutate(theoretical = bases/count)  %>%
    group_by(samplename, prep) %>% 
    do(data_frame(coverage = seq(max(df$coverage)),
                density = dpois(seq(max(df$coverage)), .$theoretical) ))%>%
    ungroup() %>%
    mutate(line_type = 'Theoretical (Poisson)') 
    
    
plot_df <- rbind(base_df, df %>%select(-baseq_count, -count)) %>%
    mutate(line_type = factor(line_type, levels = c('WGS','Theoretical (Poisson)'))) %>%
#    filter(grepl('^K12',samplename)) %>%
    tbl_df

rsqrd_df <- plot_df %>% 
    spread(line_type, density) %>% 
    set_names(c('samplename','prep','coverage','wgs','model')) %>% 
    filter(!is.na(model))  %>%
    group_by(samplename, prep) %>%
    summarize(sum_res = sum((wgs- model)^2),
              sum_var = sum((wgs - mean(wgs))^2),
              rmsd = sqrt(sum((wgs-model)^2)/length(wgs))
    ) %>%
    ungroup() %>%
    mutate(rsqrd = 1 - sum_res/ sum_var) %>%
    filter(grepl('nextera|umi',samplename)) 

rsqrd <- rsqrd_df %>%
    group_by(prep) %>%
    summarize(mean_rsqd = signif(mean(rsqrd),3),
              sd_rsqd = signif(sd(rsqrd),3)) %>%
    tbl_df

plot_df <- inner_join(plot_df, rsqrd) %>%
    mutate(prep = str_c(prep, ' (R-sqrd: ',mean_rsqd,'±',sd_rsqd,')'))
    
wgs_p <- ggplot() +
	geom_line(data = plot_df, aes(x = coverage, y = density * 100, 
	                   color = prep, size = samplename, linetype=line_type)) + 
	xlim(1,50) +
    scale_color_manual(values = c('light sky blue','salmon'))+
    scale_size_manual(guide = 'none', values = rep(1.2, length(unique(plot_df$samplename))))+
    scale_linetype_discrete(guide = guide_legend(ncol = 1))+
	theme(text = element_text(size = 25, face='bold')) +
	theme(axis.text = element_text(size = 25, face='bold')) +
	labs(x ='Level of Coverage', y = '% of Genome', 
	     color =' ', linetype = ' ') +
    theme(legend.position = c(0.7,0.5)) +
    theme(legend.key.height=unit(2,"line"))
ggsave(wgs_p, file = figurename, height=7,width=10)
message('Plotted: ', figurename)


# rsqrd_p <- ggplot(data = rsqrd_df, aes(x = prep, y = rsqrd, 
#                                        color = prep)) +
#     geom_boxplot(size=2) +
#     labs(x = ' ', y=expression(R^{2}),parse=T) +
#     theme(legend.position = 'none') +
#     theme(axis.text.y = element_text(face='bold',size = 20))+
#     theme(text = element_text(size = 18)) +
#     theme(axis.text = element_text(size = 18)) +
#     theme(axis.ticks.y = element_blank()) +
#     theme(axis.text.x = element_text(size = 18, angle = 60, hjust = 1, vjust =1)) +
#     coord_flip()
# 
# model_df <- plot_df %>% 
#     spread(line_type, density) %>% 
#     set_names(c('samplename','prep','coverage','wgs','model')) %>% 
#     filter(!is.na(model))  %>%
#     group_by(samplename, prep) %>%
#     nest() %>%
#     mutate(lm_model = map(data, ~lm(wgs~model, data = .))) %>%
#     mutate(anova_test = map(lm_model, anova)) %>%
#     mutate(summary_stat = map(anova_test, tidy))  %>%
#     unnest(summary_stat)  %>%
#     filter(term == 'Residuals') 
# 
# residual_p <- ggplot(data = model_df, aes(x = prep, y = sumsq, 
#                                        color = prep)) +
#     geom_boxplot(size=2) +
#     labs(x = ' ', y='Residual Sum of Squares') +
#     theme(legend.position = 'none') +
#     theme(axis.text.y = element_text(face='bold',size = 20))+
#     theme(text = element_text(size = 18)) +
#     theme(axis.text = element_text(size = 18)) +
#     theme(axis.ticks.y = element_blank()) +
#     theme(axis.text.x = element_text(size = 18, angle = 60, hjust = 1, vjust =1)) +
#     coord_flip()

   