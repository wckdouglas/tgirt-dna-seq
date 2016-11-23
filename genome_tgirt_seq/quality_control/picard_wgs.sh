#!/usr/bin/bash

PROJECT_PATH=${WORK}/cdw2854/ecoli_genome
BAM_PATH=$PROJECT_PATH/bamFiles/subsampled_bam
PICARD_PATH=$PROJECT_PATH/picard_results
REF_PATH=$REF/Ecoli
GENOME=$REF_PATH/k12_mg1655.fa

for BAM in $BAM_PATH/*.MarkDuplicate*.bam
do
	SAMPLE_NAME=$(basename ${BAM%.sorted.bam})
	echo  picard CollectRawWgsMetrics INPUT=${BAM} \
		OUTPUT=$PICARD_PATH/${SAMPLE_NAME}.wgs.metrics \
		REFERENCE_SEQUENCE=$GENOME \
		INCLUDE_BQ_HISTOGRAM=true
done
