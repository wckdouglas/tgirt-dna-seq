#!/usr/bin/env python

library(zoo)
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
        mutate(subsampled = ifelse(grepl('1M',filename),'Yes','No'))  %>%
        tbl_df
	message('Read ', filename)
	return(df)
}

project_path <- '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path <- str_c( project_path, '/picard_results')
figure_path <- str_c( project_path, '/figures')
figurename <- str_c( figure_path, '/wgs_plot.pdf')
table_names <- list.files(path = picard_path , pattern = '.wgs.metrics')
table_names <- table_names[grepl('1M',table_names)]
df <- table_names %>%
	map_df(read_wgs_table, picard_path) %>%
    mutate(line_type = 'WGS') %>%
    select(-subsampled) %>%
    filter(!grepl('sim',samplename)) %>%
    filter(!grepl('clustered',samplename)) %>%
    filter(samplename != 'SRR733099') %>%
    tbl_df

base_df <- df %>%
    mutate(bases = count * coverage) %>%
    group_by(samplename) %>%
    summarize(bases = sum(bases),
              count = sum(count)) %>%
    ungroup() %>%
    mutate(theoretical = bases/count)  %>%
    group_by(samplename) %>% 
    do(data_frame(coverage = seq(max(df$coverage)),
                density = dpois(seq(max(df$coverage)), .$theoretical) ))%>%
    ungroup() %>%
    mutate(line_type = 'Theoretical (Poisson)') 
    
    
plot_df <- rbind(base_df, df %>%select(-baseq_count, -count)) %>%
    mutate(line_type = factor(line_type, levels = c('WGS','Theoretical (Poisson)'))) %>%
#    filter(grepl('^K12',samplename)) %>%
    filter(!grepl('X',samplename)) %>%
    tbl_df

    
wgs_p <- ggplot() +
	geom_line(data = plot_df, aes(x = coverage, y = density, 
	                              color = samplename, linetype = line_type), 
	          size = 1.3) + 
	xlim(1,50) +
	theme(legend.position = c(0.6, 0.9)) +
#    scale_color_discrete(guide = F)+
    scale_linetype_discrete(guide = guide_legend(ncol = 1))+
	theme(text = element_text(size = 20)) +
	theme(axis.text = element_text(size = 18)) +
	labs(x ='Depth of coverage', y = '% of Genome', 
	     color =' ', linetype = ' ') 
#ggsave(wgs_p, file = figurename)
#message('Saved: ', figurename)



model_df <- base_df %>% 
    select(coverage, density,samplename) %>%
    dplyr::rename(poisson=density) %>%
    inner_join(df %>% select(density, coverage, samplename)) %>%
    filter(!grepl('[86]X',samplename)) %>%
    filter(grepl('^K12',samplename)) %>%
    group_by(samplename) %>%
    summarize(sum_res = sum((density - poisson)^2),
              sum_var = sum((density - mean(density))^2),
              rmsd = sqrt(sum((density-poisson)^2)/length(density))
              ) %>%
    ungroup() %>%
    mutate(rsqrd = 1 - sum_res/ sum_var) %>%
    mutate(enzyme = str_sub(samplename,5,6)) %>%
    mutate(enzyme = sapply(enzyme, rename_enzyme))  %>%
    group_by(enzyme) %>%
    summarize(mean_rsqrd = mean(rsqrd),
              min_rsqrd = min(rsqrd),
              max_rsqrd = max(rsqrd)) %>%
    ungroup()

library(scales)
p <- ggplot(data = model_df, aes(x = enzyme, y = mean_rsqrd, fill = enzyme)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(ymin = min_rsqrd, ymax = max_rsqrd), width = 0.4) +
    theme(legend.position = 'none') +
    labs(x = ' ', y = 'R-sqaured', title='Model fitting to Possion') +
    scale_y_continuous(limits=c(0.8,1),oob = rescale_none)
figurename <- str_c(figure_path, '/enzyme_rsqrd_plot.pdf')
ggsave(p , file = figurename)
    