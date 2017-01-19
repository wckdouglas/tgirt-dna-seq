#!/usr/bin/bash

BW_PATH=$WORK/cdw2854/plasmaDNA/genomeWPS/bigWig_files
MERGE_BW_PATH=$BW_PATH/mergedWPS
CHROM_SIZE_TSV=$REF/GRCh38/hg38_rDNA/hg38_genome.tsv

for SAMPLENAME in $(ls $BW_PATH/*.bigWig | cut -d'.' -f1 | uniq)
do
	for LENGTH in Short Long
	do
		BED_GRAPH=$MERGE_BW_PATH/$(basename ${SAMPLENAME}).bedGraph
		MERGED_BW=${BED_GRAPH/.bedGraph/.bigWig}
		echo bigWigMerge \
			${SAMPLENAME}.\*.${LENGTH}.bigWig \
			/dev/stdout \
			\| sort -k1,1 -k2,2n \
			\> $BED_GRAPH\
		\;bedGraphToBigWig \
			$BED_GRAPH \
			${CHROM_SIZE_TSV} \
			$MERGED_BW
	done
done
