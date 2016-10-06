#!/usr/bin/bash

PROJECT_PATH=${WORK}/cdw2854/ecoli_genome
PICARD_PATH=$PROJECT_PATH/picard_results
REF_PATH=$REF/Ecoli
GENOME=$REF_PATH/b_strain.fa

for BAM in $PICARD_PATH/*.MarkDuplicate.bam
do
	SAMPLE_NAME=$(basename ${BAM%.sorted.bam})
	echo  picard CollectRawWgsMetrics INPUT=${BAM} \
		OUTPUT=$PICARD_PATH/${SAMPLE_NAME}.wgs.metrics \
		REFERENCE_SEQUENCE=$GENOME \
		INCLUDE_BQ_HISTOGRAM=true
done
