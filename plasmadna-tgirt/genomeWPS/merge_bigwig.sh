#!/usr/bin/bash

BW_PATH=/scratch/02727/cdw2854/plasmaDNA/genomeWPS
MERGE_BW_PATH=$BW_PATH/mergedWPS
CHROM_SIZE_TSV=$REF/GRCh38/hg38_rDNA/hg38_genome.tsv

for SAMPLENAME in $(ls $BW_PATH/*.bigWig | cut -d'.' -f1 | uniq)
do
	for LENGTH in Short Long
	do
		echo bigWigMerge \
			${SAMPLENAME}.*.${LENGTH}.bigWig \
			/dev/stdout \
		\| bedGraphToBigWig \
			/dev/stdin \
			${CHROM_SIZE_TSV} \
			${MERGE_BW_PATH}/${SAMPLENAME}.bigWig
	done
done
