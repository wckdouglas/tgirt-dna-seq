#!/usr/bin/env Rscript

library(Hmisc)
library(parallel)
library(tidyr)
library(readr)
library(quadprog)
library(cowplot)
library(stringr)
library(dplyr)

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
		filter(coverage > 5) %>%
		mutate(periodontium = (periodontum + periodontium)/2) %>%
		select(-placenta, - breast, -embryo, -periodontum, -`spinal cord`) %>%
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

main <- function(){
	filepath <- '/scratch/02727/cdw2854/plasmaDNA/tissueContribution'
	filepath <- '/scratch/cdw2854/plasmaDNA/tissueContribution'
	figurepath <- '/scratch/02727/cdw2854/plasmaDNA/figures'
	figurepath <- '/scratch/cdw2854/plasmaDNA/figures'
	files <- list.files(path=filepath, pattern = '.tsv')
	files <- files[grepl('RNase|NT|SRR|RNA-OCD|ref',files)]

	df <- mclapply(files,readFile, filepath, mc.cores=24) %>%
#	df <- lapply(files,readFile, filepath) %>%
		do.call(rbind,.)  %>%
		tbl_df

			
	labeling <- df %>%
		group_by(samplename) %>%	
		do(data_frame(pos = cumsum(.$fraction)- .$fraction/2,
				  tissue = .$tissue,
				  fraction = .$fraction))  %>%
		filter(fraction > 0.05) %>%
		ungroup() %>%
		tbl_df

	p <- ggplot(data = df, aes(x=factor(1),fill=tissue,y=fraction)) + 
		geom_bar(stat='identity')+
		coord_polar(theta='y') + 
		geom_text(data= labeling, aes(label=tissue, y=pos, x=factor(1))) +
		theme(axis.text.x=element_blank()) +
		theme(axis.text.y=element_blank()) +
		labs(x=' ', y =' ', fill='Tissues') +
		theme(axis.line=element_blank()) +
		facet_wrap(~samplename)		

#	p <- ggplot(data = df, aes(x=samplename,fill=tissue,y=fraction)) + 
#		geom_bar(stat='identity')+
#		labs(x='Sample', y ='Fraction of contribution', fill='Tissues') +
#		theme(axis.line=element_blank()) +
#		theme(text=element_text(size=20,face='bold'))
	figurename <- str_c(figurepath,'/tissueDHS.pdf')
	ggsave(p, file=figurename)
	message('Written: ', figurename)
}

main()
