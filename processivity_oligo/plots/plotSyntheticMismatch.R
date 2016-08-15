#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(stringi)
library(tidyr)
library(cowplot)
library(plasmaDNA)

datapath <- '/scratch/cdw2854/tgirtDNA/baseMismatchCount'
figurepath <- '.'

getSampleType <- function(name){
	ifelse(!grepl('rror',name),'Raw','corrected')
}

getReplicate <- function(name){
	ifelse(!grepl('-1',name),'Replicate 1','Replicate 2')
}

readFile <- function(filename,datapath){
	if (file.info(stri_c(datapath,filename,sep='/'))$size > 0) {
		datapath %>%
		str_c(filename,sep='/') %>%
		read_tsv() %>%
		mutate(name = filename ) %>%
		mutate(refBase=toupper(refBase)) %>%
		tbl_df %>%
		return 
	}
}

files <- list.files(path = datapath,pattern='tsv')
resultDF <- lapply(files,readFile,datapath) %>%
	.[!sapply(.,is.null)] %>%
	do.call(rbind,.) %>%
	tbl_df

coverage <- resultDF %>%
	mutate(template = assignTemplate(name)) %>%
	mutate(prep = getSampleType(name)) %>% 
	mutate(enzyme = assignEnzyme(name)) %>%
	mutate(replicate = getReplicate(name)) %>%
	group_by(template, prep, enzyme, replicate) %>%
	summarize(cov = sum(coverage))  %>%
	tbl_df

df <- resultDF %>% 
	select(A,C,T,G,refBase,coverage,name) %>% 
	gather(mismatch,count,-refBase:-name) %>% 
	filter(refBase!=mismatch) %>% 
	mutate(mismatch = str_c(refBase,mismatch,sep='>')) %>%
	mutate(template = assignTemplate(name)) %>%
	mutate(prep = getSampleType(name)) %>% 
	mutate(enzyme = assignEnzyme(name)) %>%
	mutate(replicate = getReplicate(name)) %>%
	group_by(template,prep,enzyme,mismatch, replicate) %>% 
	summarize(count=sum(count)) %>%
	inner_join(coverage) %>%
	mutate(count = count/cov) %>%
	ungroup() %>%
	select(-cov) %>%
	filter(prep != 'Merge correct') %>%
	tbl_df

plotFunction <- function(replicateNum,df){
	enzyme <- unique(df$enzyme)
	coverageDF <- coverage %>% filter(replicate == replicateNum) 
    p <- ggplot(data = df[df$replicate == replicateNum,],aes(x=mismatch,y=count*1e4,fill=prep)) +
	    geom_bar(stat='identity') +
	    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
	    geom_text(data=coverageDF[coverageDF$enzyme==enzyme,],aes(x=3,y=6,label=str_c('Base-in: ',cov))) +
	    facet_grid(template~prep) +
	    panel_border() +
	    labs(y=expression(Mutation~rate~"(x"*10^-4*")"),x = ' ',fill=' ') +
	    theme(text = element_text(size=20,face='bold')) +
	    theme(legend.position='none')  
	return(p)
}

plotRep <- function(enz, df){
	df1 <- df[df$enzyme==enz,]
	ps <- lapply(unique(df1$replicate),plotFunction,df1)
	p <- plot_grid(plotlist=ps, ncol=2)
	return (p)
}

ps <- lapply(unique(df$enzyme),plotRep,df)
p <- plot_grid(plotlist=ps,labels=unique(df$enzyme),ncol=1)
figurename <- str_c(figurepath,'/mismatch.pdf')
ggsave(p,file=figurename,width=15,height=15)
message('Saved: ',figurename)


