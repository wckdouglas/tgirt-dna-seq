#!/usr/bin/env Rscript

library(cowplot)
library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(openxlsx)

table_dir <- '/stor/work/Lambowitz/cdw2854/target-seq/base_tables'
tables <- list.files(path = table_dir, pattern = '.tsv')

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
#        filter(chrom %in% c("ENST00000613854","ENST00000612966","ENST00000621411",
#                            "ENST00000616365","ENST00000359303","ENST00000356476",
#                            "ENST00000614911","ENST00000618052")) %>%
        tbl_df
    return (df)
}

transition <- c('A-to-G','G-to-A','C-to-T','T-to-C')

assignTemplate <- function(x){
    ifelse(grepl('^H',x),'RNA','DNA')
}

assignFilter <- function(x){
    ifelse(grepl('.4',x),'>4 members','>3 members')
}

assignMethod <- function(x){
    ifelse(grepl('bayesian',x),'New clustering','Traditional clustering')
}

figurename <- str_c(table_dir, '/mismatch.pdf')
merge_df <- tables %>%
    map(read_file) %>%
    reduce(rbind)  %>%
    filter(grepl('filter',samplename)) %>%
    tbl_df

clean_mismatch_df <- merge_df %>%
    select(ref_base, A,C,T,G, samplename, start) %>%
    gather(base, count, -samplename, -ref_base, -start)  %>%
    tbl_df

mismatch_df <- clean_mismatch_df %>%
    group_by(samplename, base, ref_base) %>%
    summarize(count = sum(count)) %>%
    ungroup() %>%
    mutate(mutation = str_c(ref_base, '-to-', base)) %>%
    group_by(samplename) %>%
    do(data_frame(
        count = .$count/sum(.$count),
        mutations = .$mutation,
        ref = .$ref_base,
        read = .$base
    )) %>%
    ungroup() %>%
    filter(ref!=read) %>%
    mutate(filter_data = assignFilter(samplename)) %>%
    mutate(template = assignTemplate(samplename) ) %>%
#    mutate(samplename = str_c(template, ' (',filter_data, ')')) %>%
    mutate(label_mutation = ifelse(mutations %in% transition, 'Transition','Transversion')) %>%
    tbl_df

total_base <- merge_df %>%
    group_by(samplename) %>%
    summarize(total_base = sum(coverage)) %>%
    inner_join(mismatch_df %>% select(samplename, mutations, count)) %>%
    mutate(real_count = count * total_base) %>%
    group_by(samplename) %>%
    do(data_frame(
        total_base = .$total_base,
        mutations = .$mutations,
        count = .$real_count/ sum(.$real_count)
    )) %>%
    spread(mutations,count)

error_DF <- merge_df %>%
    select(-end,-strand) %>%
    gather(base, count, -samplename, -chrom, -start, -ref_base, -coverage) %>%
    filter(ref_base != base) %>%
    group_by(samplename, start) %>%
    summarize(count = sum(count),
              coverage = unique(coverage)) %>%
    ungroup() %>%
    group_by(samplename) %>%
    summarize(count = sum(count),
              coverage = sum(coverage)) %>%
    ungroup() %>%
    mutate(error_rate = count / coverage) %>%
    tbl_df

tablename <- str_replace(figurename,'.pdf','.csv')
lambowitz_table <- merge_df %>%
    group_by(samplename) %>%
    summarize(deletions = sum(deletions),
              insertions = sum(insertions)) %>%
	mutate(indel = deletions + insertions) %>%
    inner_join(total_base)  %>%
    mutate(deletion_rate = deletions / total_base) %>%
    mutate(insertion_rate = insertions / total_base) %>%
	mutate(indel_rate = deletion_rate + insertion_rate ) %>%
    inner_join(error_DF %>% select(samplename, error_rate)) %>%
    mutate(template = assignTemplate(samplename)) %>%
    mutate(method = assignMethod(samplename)) %>%
    mutate(filter_member = assignFilter(samplename)) %>%
	mutate(mismatch_rate = error_rate - insertion_rate - deletion_rate) %>%
    gather(errors, count,-samplename) %>%
    spread(samplename, count) %>%
    mutate(errors = str_replace_all(errors,'-',' ')) %>%
    mutate(errors = as.factor(errors)) %>%
    mutate(errors = relevel(errors, 'deletions')) %>%
    mutate(errors = relevel(errors, 'insertions')) %>%
    mutate(errors = relevel(errors, 'indel')) %>%
    mutate(errors = relevel(errors, 'deletion_rate')) %>%
    mutate(errors = relevel(errors, 'insertion_rate')) %>%
    mutate(errors = relevel(errors, 'indel_rate')) %>%
    mutate(errors = relevel(errors, 'mismatch_rate')) %>%
    mutate(errors = relevel(errors, 'error_rate')) %>%
    mutate(errors = relevel(errors, 'total_base')) %>%
    mutate(errors = relevel(errors, 'template')) %>%
    mutate(errors = relevel(errors, 'filter_member')) %>%
    mutate(errors = relevel(errors, 'method')) %>%
#	mutate(errors = str_replace_all(errors, '_',' ')) %>%
    arrange(errors) %>%
    write_csv(tablename)
message('Written ',tablename)
