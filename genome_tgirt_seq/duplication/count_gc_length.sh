#!/bin/bash

BED_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id/bam_files/bed_files

for BED in $BED_PATH/*.bed
do
	OUT_FILE=${BED%.bed}.tsv
	echo "duplication\tlength\tgc_pct" > $OUT_FILE
	echo bedtools nuc -bed $BED -fi $REF/Ecoli/k12_mg1655.fa \
		\| cut -f4,5,8 \
		\| sed 1d \
		\| python parse_gc.py \
		\>\> $OUT_FILE
done
