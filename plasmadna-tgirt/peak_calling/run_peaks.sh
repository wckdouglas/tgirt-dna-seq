#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/plasmaDNA
#PROJECT_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA
BW_PATH=$PROJECT_PATH/genomeWPS/bigWig_files
BED_PATH=$PROJECT_PATH/genomeWPS/bed_files

for BW in $BW_PATH/*bigWig
do
	FILENAME=$(basename $BW)
	SAMPLENAME=$(echo $FILENAME | cut -d'.' -f1)
	CHROM=$(echo $FILENAME | cut -d'.' -f2)
	TYPE=$(echo $FILENAME | cut -d'.' -f3)
	echo python peak_calling.py \
		--in_bigwig $BW \
		--out_bed $BED_PATH/${FILENAME/.bigWig/.bed} \
		--length_type=$TYPE \
		--chrom=$CHROM
done
