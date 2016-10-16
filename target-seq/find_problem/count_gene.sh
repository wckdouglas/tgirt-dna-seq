#!/usr/bin/bash

BAM_PATH=/stor/work/Lambowitz/cdw2854/target-seq/removed_clipped

parallel samtools index {} ::: $BAM_PATH/*bam
parallel samtools idxstats {} \> {.}.tsv ::: $BAM_PATH/*bam
python make_gene_table.py
Rscript plot_genes.R
