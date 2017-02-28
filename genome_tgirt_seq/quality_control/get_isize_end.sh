#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/ecoli_genome/
BAM_PATH=$PROJECT_PATH/picard_results
OUT_PATH=$PROJECT_PATH/isize_tables
PROGRAM_PATH=/stor/home/cdw2854/ngs_qc_plot
mkdir -p $OUT_PATH

for BAM in $BAM_PATH/*.MarkDuplicate.bam 
do 
	SAMPLENAME=$(basename ${BAM%.bam})
	echo python $PROGRAM_PATH/bam_insert_size.py $BAM ${OUT_PATH}/${SAMPLENAME}.tsv
done

