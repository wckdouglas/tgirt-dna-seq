#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS
BIGWIG_PATH=$PROJECT_PATH/bigWig_files

for BW in $BIGWIG_PATH/*long.bigWig
do
	SAMPLENAME=$(basename $BW)
	CHROM=$(echo $SAMPLENAME | rev | cut -d'.' -f3 | rev)
	for GENE_REGION in gene_body tss
	do
		OUTPATH_DIR=$PROJECT_PATH/periodicity_$GENE_REGION
		for FFT in scipy R
		do
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

