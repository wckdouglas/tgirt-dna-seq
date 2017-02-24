#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/ecoli_genome/
BAM_PATH=$PROJECT_PATH/picard_results
OUTT_PATH=$PROJECT_PATH/fragment_ends
for BAM in $BAM_PATH/*.MarkDuplicate.bam 
do 
	SAMPLENAME=$(basename ${BAM%.bam})
	echo python ~/ngs_qc_plot/bam_fragment_ends.py \
		$BAM \
		$PROJECT_PATH/fragment_ends/$SAMPLENAME
done

