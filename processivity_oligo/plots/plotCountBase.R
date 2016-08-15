#!/usr/bin/env Rscript

library(cowplot)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(plasmaDNA)

getSampleType <- function(name){
	ifelse(!grepl('rror',name),'Raw','corrected')
}

datapath <- '/scratch/cdw2854/plasmaDNA/baseMismatchCount'
figurepath <- '.'
files <- list.files(path=datapath, pattern = '.count.tsv')

readFile <- function(filename, datapath){
	if (file.info(str_c(datapath,'/',filename))$size != 0){
		samplename <- str_split(filename,'\\.')[[1]][1]
		df <- datapath %>%
			str_c(filename,sep='/') %>%
			read_tsv() %>%
			gather(read,count,-ref) %>%
			mutate(name = samplename) %>%
			tbl_df
		return(df)
	}	
}

result <- lapply(files, readFile, datapath) %>%
	.[!sapply(.,is.null)] %>%
	do.call(rbind,.) %>%
	mutate(template = assignTemplate(name)) %>%
	mutate(prep = getSampleType(name)) %>% 
	mutate(enzyme = assignEnzyme(name)) %>%
	group_by(template, enzyme, prep, ref, read) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	group_by(template, enzyme, prep) %>%
	do(data_frame(coverage = sum(.$count),
				  count = .$count,
				  ref = .$ref,
				  read = .$read)) %>%
	ungroup() %>%
	mutate(mismatch = str_c(ref,read,sep='>')) %>%
	mutate(count = count / coverage) %>%
	filter(ref!=read) %>%
	tbl_df 

	
plotFunction<-function(enzyme,df){
    ggplot(data=df[df$enzyme==enzyme,],aes(x=mismatch,y=count*1e4,fill=prep)) +
	    geom_bar(stat='identity') +
	    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
	    facet_grid(template~prep) +
	    panel_border() +
	    labs(y=expression(Mutation~rate~"(x"*10^-4*")"),x = ' ',fill=' ') +
	    theme(text = element_text(size=20,face='bold')) +
	    theme(legend.position='none')  
}
ps <- lapply(unique(result$enzyme),plotFunction,result)
p <- plot_grid(plotlist=ps,labels=unique(result$enzyme),ncol=1)
figurename <- str_c(figurepath,'/mismatchCount.pdf')
ggsave(p,file=figurename,width=15,height=15)
message('Saved: ',figurename)


