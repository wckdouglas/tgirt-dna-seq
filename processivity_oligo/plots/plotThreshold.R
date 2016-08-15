#!/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)

p <- read_csv('try.dat',col_names=c('loglik','cov','match'),col_type='nnn') %>%
    filter(match!=0) %>%
    mutate(ratio = match/cov) %>%
    filter(loglik < 10) %>%
    group_by(loglik, ratio,cov, match) %>%
    summarize(count = n()) %>%
    ggplot(aes(y= loglik, x=ratio)) +
        geom_point(alpha=0.7,aes(color=cov,size=count)) +
        scale_x_continuous(breaks=seq(0,1,0.05)) +
        scale_y_continuous(breaks=seq(0,5,0.5)) +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
        geom_smooth(method='lm') +
        geom_hline(y=4, color='blue') +
        scale_color_gradient(high='red',low='yellow') +
        labs(x='Correct base ratio',y = 'Log likelihood ratio', color='Base coverage', size = 'Number of occurance')
ggsave(p,file='try.pdf')

