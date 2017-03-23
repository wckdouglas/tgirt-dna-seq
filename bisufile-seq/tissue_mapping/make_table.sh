#!/bin/bash

REF_TABLE=$REF/hg19/methylation/methyl_table.bed
SAMPLE_TABLE=$WORK/cdw2854/bisufite_seq/methyl_calling/P13B_mix_S1_CpG.meth.bedGraph
OUT_PATH=$WORK/cdw2854/bisufite_seq/tissue_table
OUT_TABLE=$OUT_PATH/$(basename ${SAMPLE_TABLE%.bedGraph}.bedGraph)
mkdir -p $OUT_PATH

echo bedtools intersect -b $REF_TABLE -a $SAMPLE_TABLE -wb \
	\| cut -f1-4,8 \
	\| datamash groupby 5 mean 4 \
	\| awk \''{print $1,$2*100}'\' OFS=\'\\t\' \
	\| tr \'':-'\' \''\t\t'\' \
	\> $OUT_TABLE
