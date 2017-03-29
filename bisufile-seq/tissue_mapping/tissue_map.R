#!/usr/bin/env Rscript

library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(cowplot)
library(quadprog)
library(forcats)


solveQP <- function(Xmat, Y){
    #http://stats.stackexchange.com/questions/21565/how-do-i-fit-a-constrained-regression-in-r-so-that-coefficients-total-1
    Rinv <- solve(chol(t(Xmat) %*% Xmat))
    C <- cbind(rep(1,ncol(Xmat)), diag(ncol(Xmat)))
    b <- c(1,rep(0,ncol(Xmat)))
    d <- t(Y) %*% Xmat  
    solution <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    return(solution$solution)
}

tissue_ref_table <- '/stor/work/Lambowitz/ref/hg19/methylation/methyl_tissue_table.tsv' %>% 
    read_tsv()%>%
    set_names(make.names(names(.))) %>%
    tbl_df

df <- '/stor/work/Lambowitz/cdw2854/bisufite_seq/tissue_table/P13B_mix_S1_CpG.bedGraph' %>%
    read_tsv(col_names = c('chrom','start','end','methyl_density')) %>%
    mutate(Genomic.location = str_c(chrom,':',start,'-',end)) %>%
    inner_join(tissue_ref_table, by = c("Genomic.location"))  %>%
    select(-chrom,-start,-end,-Placenta,-Mean)

dm <- select(df,1,3:15) %>%
    data.matrix()

result <- solveQP(dm[,-1], dm[,1])
qdf <- data_frame(tissues = colnames(dm[,-1]),
           fraction =result) %>%
    mutate(tissue_type = case_when(
            grepl('B.cells|T.cells',.$tissues) ~ "Lymphocytes",
            TRUE~str_replace(.$tissues,'\\.',' ')
    )) %>%
    group_by(tissue_type) %>%
    summarize(fraction = sum(fraction)) %>%
    mutate(tissue_type_annot = str_c(tissue_type,': ', signif(fraction*100,2),'%')) %>%
    mutate(tissue_type_annot = fct_reorder(tissue_type_annot, -fraction)) %>%
    mutate(tissue_type = fct_reorder(tissue_type, -fraction)) %>%
    arrange(tissue_type_annot) %>%
    mutate(pos = 1 - cumsum(fraction) + fraction/2)  

colors <- c('lightgoldenrod2','khaki1','purple4','orange','brown','mediumseagreen','pink',
            'darkgreen','slateblue','red','skyblue','darkblue')
pie_methyl <- ggplot(data=qdf, aes(x = factor(1),y=fraction, fill=tissue_type_annot))+
    geom_bar(stat='identity')+
    geom_text(data = qdf %>% filter(fraction > 0.1), aes(x= factor(1), y=pos, label = tissue_type), size=5) +
    theme(axis.text.x = element_text(angle=0)) +
    coord_polar(theta = "y")  +
    scale_fill_manual(values=colors) +
    labs(x = ' ',y = ' ', fill=' ')+
    theme_void() +
    theme(legend.text = element_text(size=13, family='Arial')) +
    theme(legend.key.height  = unit(1,'line'))  +
    theme(legend.position = c(1.3,0.5))
p<-ggdraw() +
    draw_plot(pie_methyl, 0, 0, 0.6, 0.6)+
    draw_plot_label(letters[1:2], c(0,0),c(1,0.45), size=25)
figurepath <- '/stor/work/Lambowitz/cdw2854/bisufite_seq/figures'
figurename <- str_c(figurepath, '/methyl_pie.pdf')
ggsave(p, file = figurename, width = 7, height=12)
message('Plotted: ', figurename)

    
