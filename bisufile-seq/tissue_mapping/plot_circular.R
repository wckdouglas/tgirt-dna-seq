#!/usr/bin/env Rscript

library(rtracklayer)
library(ggbio)
library(readr)
library(dplyr)
library(cowplot)
library("extrafont")
loadfonts() 

#make genome chromosome ranges
source('/stor/home/cdw2854/tgirt-dna-seq/bisufile-seq/tissue_mapping/tissue_map.R')
genome <- '/stor/work/Lambowitz/ref/hg19/Sequence/WholeGenomeFasta/hg19.genome' %>%
    read_tsv(col_names = c('seqname','seqlength')) %>%
    filter(!seqname %in% c('chrM','chrX','chrY')) %>%
    mutate(number = str_extract(seqname,'[0-9]+'))%>%
    arrange(as.numeric(number))
hg19_genome <- GRanges(seqnames=genome$seqname, ranges=IRanges(start=1,end=genome$seqlength))
seqinfo(hg19_genome)@seqlengths <- genome$seqlength
genome(hg19_genome) <- 'hg19'


#read annotated data binned to 500bp and extracted from markers 
bed_graph_methyl <- '/stor/work/Lambowitz/cdw2854/bisufite_seq/tissue_table/P13B_mix_S1_CpG.meth.bedGraph'
mbg<-import(bed_graph_methyl, format='bedGraph', genome='hg19')
seqlevels(mbg) <- genome$seqname


# read methylation sites
methyl_marker_table <- '/stor/work/Lambowitz/ref/hg19/methylation/methyl_table.bed'
tl<-qdf$tissue_type 
marker_II <- 'Type II biomarkers'
tl <- c(as.vector(tl),marker_II)
marker_annotation <- read_tsv(methyl_marker_table,
                              col_names=c('chrom','start','end','genome_coord','tissues'))  %>%
    #filter(tissues!='none') %>%
    mutate(tissue_type = case_when(
            grepl('B.cells|T.cells',.$tissues) ~ "Lymphocytes",
            grepl('Lung',.$tissues) ~'Lungs',
            .$tissues=='none' ~marker_II,
            TRUE~str_replace(.$tissues,'\\.',' ')
        )) %>%
    mutate(tissue_alpha = case_when(
            .$tissues=='none' ~0.1,
            TRUE~1
        )) %>%
    mutate(tissue_type = factor(tissue_type, level = tl)) %>%
    arrange(rev(tissue_type)) %>%
    tbl_df
methyl_marker <- GRanges(marker_annotation$genome_coord)
mcols(methyl_marker) <- marker_annotation %>% 
    select(tissues, tissue_type, tissue_alpha)
seqlevels(methyl_marker) <- genome$seqname
seqinfo(methyl_marker) <- seqinfo(hg19_genome)


p <- ggbio() + 
    circle(methyl_marker[mcols(methyl_marker)$tissue_type=='Type II biomarkers'], 
           geom='rect', color='grey72',fill='grey72') + # add tissue annotation
    circle( methyl_marker[mcols(methyl_marker)$tissue_type!='Type II biomarkers'], 
            geom='rect', 
            aes(color=tissue_type,fill=tissue_type, alpha=tissue_alpha)) + # add tissue annotation
    scale_alpha_continuous(guide=F) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    circle(mbg, geom='bar', aes(y = score), color = 'rosybrown') + # add methylation density bar
    circle(hg19_genome, geom = "ideo", fill = "slategrey") + # add chromosome box
    circle(hg19_genome, geom = "scale", size = 6) + # add base pair length scale
    circle(hg19_genome, geom = "text", aes(label = seqnames), vjust = -1, size = 10) + # add chrome name 
    labs(color=' ', fill= ' ') +
    theme(legend.text = element_text(size=25, family='Arial')) +
    theme(legend.key.height = unit(1.5,'line')) +
    theme(legend.position = c(0.5,0.5))
figurepath <- '/stor/work/Lambowitz/cdw2854/bisufite_seq/figures'
figurename <- str_c(figurepath, '/methyl_genome_cicular.pdf')
pdf(figurename, height=13, width=13)
p
dev.off()
message('Plotted: ', figurename)