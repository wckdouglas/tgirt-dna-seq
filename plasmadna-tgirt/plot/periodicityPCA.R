#!/bin/env Rscript

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(cowplot)
library(RNAngs)

datapath <- '/Users/wckdouglas/plasmaDNA/results/plasmaResults/periodicTable'

readTable <- function(filename, datapath){
	datapath %>%
	str_c(filename, sep='/') %>%
	read_tsv()  %>%
#	mutate(periodicity = ifelse(periodicity == -1, 0, periodicity)) %>%
	filter(periodicity!=-1) %>%
	data.table::setnames('periodicity',str_replace(filename,'.tsv','')) %>%
	return
}

files = list.files(path=datapath, pattern= 'tsv')
df <- lapply(files, readTable,datapath) %>%
	Reduce(full_join,.)  
df[is.na(df)] = 0
df <- df %>%
	mutate(rowsum = apply(.[,-1:-3], 1, sum)) %>%
	mutate(type = groupGeneType(geneType)) %>%
	filter(rowsum > 0) %>%
	select(-rowsum)  %>%
	tbl_df

plotFunction <- function(dataname, dataTable){ 
	ggplot(data=subset(dataTable,sample == dataname), aes(x=periodicity, fill=type, color=type)) + 
		scale_x_log10(breaks=seq(50,1000,50))  +
		geom_density(alpha=0.5) +
		theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

}
dataTable <- df %>% 
	gather(sample, periodicity, -geneID:-geneType, -type) %>%
	filter(type %in% c('protein_coding','tRNA','snoRNA','lincRNA'))
samplelist <- c('NT','PD2','RNA-OCD','SRR2130051')
ps <- lapply(samplelist,plotFunction, dataTable) 
p<-plot_grid(plotlist=ps, labels=samplelist)
ggsave(p,file='explore.pdf', width=18.1, height=10.2)


df <- df %>% 
	select(-SRR2130052,-refA1) %>% 
	tbl_df
pca <- prcomp(df[,c(-1:-3,-ncol(df))], center=T) 
pcaDF <- pca$x %>%
	data.frame(df[,c(1:3,ncol(df))],.) %>%
	mutate(tRNA = ifelse(type=='tRNA','yes','no')) %>%
	filter(! type %in% c('Mt','IG','TR')) %>%
	tbl_df

p <- ggplot(data=pcaDF, aes(x=PC1, y = PC2, color = type)) +
	geom_point(alpha=0.5)
ggsave(p, file='pca.pdf')

kmnDF <- df %>% 
	mutate(center = kmeans(.[,-1:-3],2)$cluster)  %>%
	tbl_df 
p <- ggplot(aes(x=NT,y=SRR2130051, shape=factor(center), color=type), data=kmnDF) +
	geom_point() +
	scale_x_log10() +
	scale_y_log10()
