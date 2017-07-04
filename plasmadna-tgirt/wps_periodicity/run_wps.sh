#!/bin/bash

PROJECT_PATH=$SCRATCH/plasmaDNA/genomeWPS
BIGWIG_PATH=$PROJECT_PATH/bigWig_files

for BW in $(ls $BIGWIG_PATH/*Long.bigWig | grep '^P' | grep -v 'P1022_1113_13_1016_mix_unique\|SRR2130051')
do
	SAMPLENAME=$(basename $BW)
	CHROM=$(echo $SAMPLENAME | rev | cut -d'.' -f3 | rev)
	for GENE_REGION in gene_body tss
	do
		OUTPATH_DIR=$PROJECT_PATH/periodicity_$GENE_REGION
		for FFT in scipy R
		do
			SAMPLENAME=$(echo $SAMPLENAME|cut -d'.' -f1)
			OUTPATH=$OUTPATH_DIR/$FFT
			mkdir -p $OUTPATH
			echo python wps_TSS.py \
				--in_bigwig=$BW \
				--chrom=$CHROM \
				--out_bed=$OUTPATH/${SAMPLENAME}.${CHROM}.bed \
				--fft_type=$FFT \
				--gene_bed=$REF/GRCh38/Bed_for_counts_only/protein.bed \
				--gene_regions=$GENE_REGION
		done
	done
done

