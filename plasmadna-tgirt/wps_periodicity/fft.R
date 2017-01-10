d <- read.table('/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/a',
                as.is=T,header=F)                                                                                                                                                                       
                                                                                                                                                                                                                                             
tseries <- as.ts(d$V3)                                                                                                                                                                                                                
fseries <- filter(c(window(tseries,1,300),tseries),
                  filter=1/seq(5,100,4),
                  method="recursive")[-(1:300)]                                                                                                                                     
demean_series <- fseries-mean(fseries,trim=0.1)                                                                                                                                                                                                  
                                                                                                                                                                                                                                             
resA1 <- spec.pgram(demean_series,pad=0.3,tap=0.3,
                    span=2,plot=F,detrend=TRUE,
                    demean=TRUE)                                                                                                                                                        
resA1$freq <- 1/resA1$freq           
