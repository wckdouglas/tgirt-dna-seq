#!/usr/bin/env Rscript

library(biomaRt)
library(dplyr)
library(readr)
library(stringr)

bedfile <- '/corral-repl/utexas/2013lambowitz/Ref/hg19/hisat2/grch37_snp/genes.bed'
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
				host="grch37.ensembl.org", 
				path="/biomart/martservice", 
				dataset="hsapiens_gene_ensembl")
attributes <- c('chromosome_name','start_position','end_position',
				'ensembl_gene_id','status','strand','gene_biotype','external_gene_name')

results<- getBM(attributes = attributes, mart = mart) %>%
	tbl_df() %>%
	filter(chromosome_name %in%  c(1:22, 'X', 'Y', 'MT')) %>%
	mutate(status=0 ) %>%
	mutate(strand = ifelse(strand >0, '+','-')) %>%
#	mutate(chromosome_name = ifelse(chromosome_name=='MT','M',chromosome_name)) %>%
#	mutate(chromosome_name = str_c('chr',chromosome_name)) %>%
	write_tsv(bedfile,  col_names=F)
