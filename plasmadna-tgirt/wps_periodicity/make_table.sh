#!/bin/bash

DATAPATH=/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/periodicity_tables
BEDPATH=$DATAPATH/bed_files
mkdir -p $BEDPATH

for TSV in $DATAPATH/*.tsv
do
	SAMPLE=$(basename ${TSV%.tsv})
	cat $TSV \
		| sed '1 d' \
		| bedtools intersect -a - \
			-b $REF/GRCh38/Bed_for_counts_only/protein.bed \
			-wb  \
		| awk '{print $8,$9,$10,$11,$12,$13,$14,$15,$5,$6}' OFS='\t' \
		| datamash -g 1,2,3,4,5,6,7,8 median 9,10 \
		> $BEDPATH/${SAMPLE}.bed
done
