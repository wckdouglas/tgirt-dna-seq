#!/usr/bin/env Rscript

library(rtracklayer)
library(ggbio)
library(readr)

#make genome chromosome ranges
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
tl<-qdf$tissue_type %>% levels
marker_annotation <- read_tsv(methyl_marker_table,
                              col_names=c('chrom','start','end','genome_coord','tissues'))  %>%
    mutate(tissue_type = case_when(
            grepl('B.cells|T.cells',.$tissues) ~ "Lymphocytes",
            grepl('Lung',.$tissues) ~'Lungs',
            TRUE~str_replace(.$tissues,'\\.',' ')
        )) %>%
    mutate(tissue_type = factor(tissue_type, level = tl))
    tbl_df
methyl_marker <- GRanges(marker_annotation$genome_coord)
mcols(methyl_marker) <- marker_annotation %>% select(tissues, tissue_type)
seqlevels(methyl_marker) <- genome$seqname
seqinfo(methyl_marker) <- seqinfo(hg19_genome)


colors <- RColorBrewer::brewer.pal(12,'Paired')
colors <- c(colors, 'black')
p <- ggbio() + 
    circle(methyl_marker, geom='rect', aes(color=tissue_type,fill=tissue_type)) + # add tissue annotation
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    circle(mbg, geom='bar', aes(y = score), color = 'rosybrown') + # add methylation density bar
    circle(hg19_genome, geom = "ideo", fill = "slategrey") + # add chromosome box
    circle(hg19_genome, geom = "scale", size = 6) + # add base pair length scale
    circle(hg19_genome, geom = "text", aes(label = seqnames), vjust = -1, size = 10) + # add chrome name 
    labs(color=' ', fill= ' ') +
    theme(legend.text = element_text(size=25)) +
    theme(legend.key.height = unit(1.5,'line')) +
    theme(legend.position = c(0.5,0.5))
figurepath <- '/stor/work/Lambowitz/cdw2854/bisufite_seq/figures'
figurename <- str_c(figurepath, '/methyl_genome_cicular.pdf')
pdf(figurename, height=15, width=15)
p
dev.off()
message('Plotted: ', figurename)