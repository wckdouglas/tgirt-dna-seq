#!/usr/bin/env Rscript

library(cowplot)
library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)


read_file <- function(tablename){
    samplename <- str_replace(tablename, '.tsv','')
    df <- table_dir %>%
        str_c(tablename, sep='/') %>%
        read_tsv() %>%
        mutate(samplename = samplename) %>%
        filter(strand == '+') %>%
        filter(start > 174) %>%
        filter(start < 370) %>%
        filter(chrom == 'ENST00000621411')  %>%
     #   filter(chrom %in% c("ENST00000621411",
     #                       "ENST00000616365","ENST00000356476")) %>%
        tbl_df
    return (df)
}

transition <- c('A -> G','G -> A','C -> T','T -> C')

assignTemplate <- function(x){
    ifelse(grepl('[0-9]+R|[WH]R_|_R_',x),'RNA','DNA')
}

plot_mismatch <- function(cleaned_mismatch_df, base_count){
    mismatch_df <- cleaned_mismatch_df %>%
        dplyr::rename(read = base) %>%
        group_by(samplename, read, ref_base) %>%
        summarize(count = sum(count)) %>%
        ungroup() %>% 
        mutate(mutations = str_c(ref_base, ' -> ', read)) %>%
        inner_join(base_count, by = c('ref_base', 'samplename')) %>%
        mutate(base_fraction = count / base_total) %>%
        filter(ref_base!=read) %>%
#        do(data_frame(
#            average_mut = mean(.$base_fraction),
#            minimum_mut = min(.$base_fraction),
#            maximum_mut = max(.$base_fraction)
#        )) %>%
        mutate(template = assignTemplate(samplename)) %>%
        mutate(label_mutation = ifelse(mutations %in% transition, 'Transition','Transversion')) %>%
        mutate(condition = case_when(grepl('^SS',.$samplename) ~ 'SuperScript II',
                                     grepl('^Pl|^HR|^420',.$samplename) ~ '420mM',
                                     grepl('^200',.$samplename) ~ '200mM',
                                     grepl('^W[DR]',.$samplename) ~ '75mM KCl')) %>%
        group_by(template, condition, label_mutation, mutations,ref_base,read) %>%
        summarize(base_fraction = mean(base_fraction)) %>%
        ungroup() %>%
        tbl_df

#    mismatch_p <- ggplot(data = mismatch_df, aes(x = mutations, y = average_mut, fill = label_mutation)) +
    mismatch_p <- ggplot(data = mismatch_df, aes(x = mutations, y = base_fraction, fill = label_mutation)) +
        geom_bar(stat='identity') +
#        geom_errorbar(aes(ymin = minimum_mut, ymax = maximum_mut), width=0.25) +
#        facet_grid(template~ref, scale='free') +
        facet_grid(template~condition, scale='free') +
        labs(x = ' ', y = 'Fractions', fill = ' ') +
        theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust =0.5)) +
        theme(text = element_text(size = 20, face='bold')) +
        theme(axis.text = element_text(size = 18, face='bold')) +
        panel_border() +
        theme(legend.position = 'bottom')
	return (mismatch_p)
}


plot_pos_mismatch <- function(clean_mismatch_df){
    pos_mismatch <- clean_mismatch_df %>%
		group_by(samplename, start, ref_base) %>%
		do(data_frame(base = .$base,
					  percentage = .$count/sum(.$count))) %>%
		ungroup() %>%
        filter(ref_base != base) %>%
        group_by(samplename, start, ref_base, base) %>%
        summarize(percentage = sum(percentage) * 100 )   %>%
        ungroup() %>%
        mutate(template = assignTemplate(samplename) ) %>%
        group_by(template, start, ref_base, base) %>%
        do(data_frame(
            average_mut = mean(.$percentage),
            minimum_mut = min(.$percentage),
            maximum_mut = max(.$percentage)
        )) %>%
        ungroup() %>%
        mutate(mutations = str_c(ref_base, base, sep='->')) %>%
        mutate(gene_name = 'HIST1H3B') %>%
        tbl_df
    
    library(RColorBrewer)
    n <- 60
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    pos_p <- ggplot(data = pos_mismatch, aes(x = start, y = average_mut, fill=mutations)) +
        geom_bar(stat='identity') +
        facet_grid(template~., scale = 'free_y') +
        theme(text = element_text(size = 20, face='bold')) +
        theme(axis.text = element_text(size = 18, face='bold')) +
        labs(x = 'Position on GOI', y = '% of Mismatch', fill= ' ') +
        theme(legend.position = 'top') + 
        guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
        scale_fill_manual(values = col_vector) +
        panel_border()
	return(pos_p)
}

table_dir <- '/stor/work/Lambowitz/cdw2854/target-seq/base_tables'
tables <- list.files(path = table_dir, pattern = '.4.tsv')
figurename <- str_c(table_dir, '/mismatch_target.pdf')

merge_df <- tables %>%
	map(read_file) %>%
	reduce(rbind)   %>%
    filter(grepl('_u_',samplename)) %>%
    filter(!grepl('200',samplename)) %>%
    tbl_df

base_count <-  merge_df %>% 
    group_by(samplename, ref_base) %>% 
    summarize(base_total = sum(coverage))
	
cleaned_mismatch_df <- merge_df %>%
    select(ref_base, A,C,T,G, samplename, start) %>%
	gather(base, count, -samplename, -ref_base, -start)  %>%
	tbl_df

mismatch_p <- plot_mismatch(cleaned_mismatch_df , base_count)
pos_p <- plot_pos_mismatch(cleaned_mismatch_df)
    
#p <- ggdraw() +
#    draw_plot(cov_p, 0, 0.666, 1, 0.333) +
#    draw_plot(pos_p, 0.02, 0.333, 0.98, 0.333) +
#    draw_plot(mismatch_p, 0, 0, 1, 0.333) 
#    draw_plot_label(label = c('(B)','(C)','(D)'), 
#                    x = c(0,0,0),
#                    y = c(1,0.666,0.333), 
#                    size = 20)
p <- plot_grid(pos_p, mismatch_p, ncol=1)
ggsave(p , file = figurename, width = 10, height = 15)
message('Plotted: ',figurename)


