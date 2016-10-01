#!/usr/bin/bash

PROJECT_PATH=${WORK}/cdw2854/genomeDNA
PICARD_PATH=$PROJECT_PATH/picard_results
REF_PATH=$REF/hg19/Sequence/WholeGenomeFasta
GENOME=$REF_PATH/genome.fa

for BAM in $PICARD_PATH/*.sorted.bam
do
	SAMPLE_NAME=$(basename ${BAM%.sorted.bam})
	echo  picard CollectRawWgsMetrics INPUT=${BAM} \
		OUTPUT=$PICARD_PATH/${SAMPLE_NAME}.wgs.metrics \
		REFERENCE_SEQUENCE=$GENOME \
		INCLUDE_BQ_HISTOGRAM=true
done
