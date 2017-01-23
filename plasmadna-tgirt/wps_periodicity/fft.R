library(dplyr)

d <- read.table('/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/try.wig',
                as.is=T,header=F)                                                                                                                                                                       
                                                                                                                                                                                                                                             
tseries <- as.ts(d$V2)                                                                                                                                                                                                                
fseries <- stats::filter(c(window(tseries,1,300),tseries),
                  filter=1/seq(5,100,4),
                  method="recursive")[-(1:300)]                                                                                                                                     
demean_series <- fseries-mean(fseries,trim=0.1) 


plot(fseries,type='l',col='green');
lines(demean_series, col='blue');
lines(tseries,col='red')                                                                                                                                                                                


resA1 <- spec.pgram(demean_series,pad=0.3,tap=0.3,
                    span=2,plot=F,detrend=TRUE,
                    demean=TRUE)                                                                                                                                                        
resA1$freq <- 1/resA1$freq           

freq_df <- data.frame(freq = resA1$freq, spec = resA1$spec) %>%
    filter(freq <= 280, freq >= 120) %>%
    tbl_df
library(ggplot2)
p<-ggplot(data=freq_df, aes(x = freq, y = spec)) + 
    geom_line() +
    scale_x_continuous(name = seq(120,280,10),breaks=seq(120,280,10))+
    theme(axis.text.x = element_text(angle=90))