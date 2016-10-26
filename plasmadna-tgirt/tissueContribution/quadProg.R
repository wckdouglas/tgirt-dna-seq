#!/usr/bin/env Rscript

library(Hmisc)
library(purrr)
library(tidyr)
library(readr)
library(quadprog)
library(cowplot)
library(stringr)
library(dplyr)
library(RColorBrewer)

solveQP <- function(Xmat, Y){
	#http://stats.stackexchange.com/questions/21565/how-do-i-fit-a-constrained-regression-in-r-so-that-coefficients-total-1
	Rinv <- solve(chol(t(Xmat) %*% Xmat))
	C <- cbind(rep(1,ncol(Xmat)), diag(ncol(Xmat)))
	b <- c(1,rep(0,ncol(Xmat)))
	d <- t(Y) %*% Xmat  
	solution <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
	return(solution$solution)
}

readFile <- function(filename, filepath){
	message('reading ', filename)
	df <- filepath %>%
		str_c(filename,sep='/') %>%
		read_tsv()  %>%
		filter(coverage > 0) %>%
		mutate(periodontium = (periodontum + periodontium)/2) %>%
		select(-placenta, - breast, -embryo, -periodontum,-periodontium, -`spinal cord`,-esophagus) %>%
		tbl_df

	dataMat <- df[,-c(1:6,ncol(df))]
	Xmat <- as.matrix(dataMat)
	Y <- df[,ncol(df)]
	solution <- solveQP(Xmat, Y/sum(Y))
	message('Solved QP')
	solDT <- data.frame(fraction = solution, tissue = names(dataMat)) %>%
		mutate(tissue = capitalize(as.character(tissue))) %>%
		filter(fraction > 0 ) %>%
		mutate(samplename = str_split(filename,'\\.')[[1]][1]) %>%
		mutate(fraction = fraction/sum(fraction)) %>%
		mutate(samplename = str_replace(filename,'.tsv','')) %>%
		tbl_df
	return(solDT)
}

filepath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/tissueContribution'
figurepath <- '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures'
files <- list.files(path=filepath, pattern = '.tsv')

df <- files %>%
    map(readFile, filepath) %>%
	reduce(rbind)  %>%
    mutate(samplename = ifelse(grepl('SRR',samplename),'ssDNA-seq\n(ref 1)','TGIRT-seq')) %>%
    filter(fraction > 0.01) %>% 
    group_by(samplename) %>% 
    do(data_frame(tissue = .$tissue,
                  fraction = .$fraction/sum(.$fraction))) %>%
	tbl_df

palette <-  brewer.pal(10,'Paired')[c(6,5,4,3,2,1,7,8,9,10)]
p <- ggplot(data = df, aes(x=samplename,fill=tissue,y=fraction)) + 
	geom_bar(stat='identity')+
	labs(y =' ', x=' ',fill=' ') +
	theme(axis.line=element_blank()) +
	theme(text=element_text(size=20,face='bold')) +
	theme(axis.text.x=element_text(size=20,face='bold')) +
	theme(axis.text.y=element_text(size=20,face='bold')) +
    coord_flip() +
    theme(legend.position='bottom') +
    scale_fill_manual(values = palette ) +
    panel_border()
figurename <- str_c(figurepath,'/tissueDHS.pdf')
ggsave(p, file=figurename, width = 10, height=4)
message('Written: ', figurename )
