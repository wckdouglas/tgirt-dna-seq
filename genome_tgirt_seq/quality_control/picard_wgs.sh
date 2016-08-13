#!/usr/bin/bash

PROJECT_PATH=/stor/scratch/Lambowitz/cdw2854/plasmaDNA
BAM_PATH=$PROJECT_PATH/bamFiles
PICARD_PATH=$PROJECT_PATH/picard_results
REF_PATH=$REF/GRCh38/hg38_rDNA

for BAM in $BAM_PATH/*bam
do
	SAMPLE_NAME=$(basename ${BAM%.bam})
	echo picard  SortSam INPUT=$BAM \
		OUTPUT=/dev/stdout \
		SORT_ORDER=coordinate \
	  \| picard CollectRawWgsMetrics INPUT=/dev/stdin \
		OUTPUT=$PICARD_PATH/${SAMPLE_NAME}.wgs.metrics \
		REFERENCE_SEQUENCE=$REF_PATH/genome_rDNA.fa \
		INCLUDE_BQ_HISTOGRAM=true
done
