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
    samplename <- str_replace_all(filename,'.subsampled.wgs_metrics','')
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
table_names <- list.files(path = picard_path , pattern = '.subsampled.wgs.metrics')
    
df <- table_names %>%
	map(read_wgs_table, picard_path) %>%
    reduce(rbind) %>%
    mutate(line_type = 'WGS') %>%
#    select(-subsampled) %>%
    filter(grepl('^75|UMI', samplename)) %>%
#    filter(grepl('nextera|UMI|NEB|kh|kq', samplename)) %>%
    filter(grepl('nextera|UMI',samplename))  %>%
    filter(grepl('nextera|umi2',samplename)) %>%
    mutate(prep = case_when(grepl('nextera',.$samplename) ~ 'Nextera~XT',
                            grepl('pb',.$samplename) ~ 'Pacbio',
                            grepl('sim',.$samplename) ~ 'Covaris Sim',
                            grepl('SRR',.$samplename) ~ 'Covaris SRR',
                            grepl('UMI',.$samplename) ~ 'TGIRT-seq 13N direct ligation',
                            grepl('kh|kq',.$samplename) ~ 'TGIRT-seq Covaris',
                            grepl('NEB',.$samplename) ~ 'TGIRT-seq Fragmentase')) %>%
    filter(!is.na(prep)) %>%
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
    mutate(rsqrd = 1 - sum_res/ sum_var) 

rsqrd <- rsqrd_df %>%
    group_by(prep) %>%
    summarize(mean_rsqd = signif(mean(rsqrd),3),
              sd_rsqd = signif(sd(rsqrd),3)) %>%
    tbl_df

plot_df <- inner_join(plot_df, rsqrd) %>%
    filter(!grepl('clust', prep)) %>%
    mutate(prep = ifelse(grepl('13N', prep),'TGIRT-seq', prep)) %>%
    mutate(clustered = ifelse(grepl('cluster',samplename),' Clustered','')) %>%
    mutate(prep = str_c(prep, clustered)) %>%
    mutate(prep = str_c(prep, '~(R^{2}:', signif(mean_rsqd,3),'%+-%',signif(sd_rsqd,1),')')) %>%
    tbl_df
preps <- plot_df$prep %>% unique

    
colors <- c('salmon','black')
wgs_p <- ggplot() +
	geom_line(data = plot_df %>% filter(line_type == 'WGS'), 
	          aes(x = coverage, y = density * 100, 
	              color = prep, group=samplename), alpha = 0.5, linetype=1) + 
    geom_line(data = plot_df %>% filter(line_type != 'WGS'), 
	          aes(x = coverage, y = density * 100, 
	              color = prep, group=samplename), alpha = 0.5, linetype=2) +
	xlim(1,30) +
    scale_color_manual(values = colors) +
    scale_linetype_discrete(guide = guide_legend(ncol = 1))+
    theme(text = element_text(size = 25, face='bold')) +
	theme(axis.text = element_text(size = 25, face='bold')) +
	labs(x ='Level of Coverage', y = '% of Genome', color =' ') +
    theme(legend.position = 'none') +
    annotate(geom='text',x=10,y=12,label=preps[1],parse=T, 
           hjust = 0, color = colors[1], size = 8) +
    annotate(geom='text',x=10,y=11,label=preps[2],parse=T, 
           hjust = 0, color = colors[2], size = 8) +
    geom_segment(aes(x = 18, xend = 20, y = 8, yend = 8), linetype=1,size = 2) +
    geom_segment(aes(x = 18, xend = 20, y = 7, yend = 7), linetype=2,size = 2) +
    annotate(geom='text', x = 21, y = 8, label = 'Experimental', size = 8, hjust =0) +
    annotate(geom='text', x = 21, y = 7, label = 'Theoretical (Poisson)', size = 8, hjust=0)
ggsave(wgs_p, file = figurename, height=8,width=9)
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

   