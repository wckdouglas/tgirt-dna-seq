#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome
BAM_PATH=$PROJECT_PATH/picard_results
MISMATCH_PATH=$PROJECT_PATH/mismatch_profiles
mkdir -p $MISMATCH_PATH

for BAM in $BAM_PATH/K12_UMI*_umi2id.MarkDuplicate.*.subsampled.bam \
			$BAM_PATH/75*nextera*MarkDuplicate.*.subsampled.bam \
			$BAM_PATH/K12_UMI*clustered*MarkDuplicate.*.subsampled.bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo mismatch_profile.py \
			-i $BAM \
			-l 75 \
			-o $MISMATCH_PATH/$SAMPLENAME \
			-q 30 \
			-n 10000
	echo samtools mpileup -f $REF/Ecoli/k12_mg1655.fa  $BAM \
		\| pileup_to_bed.py -i - -q 30 -c 1 \
		\> $MISMATCH_PATH/$SAMPLENAME.tsv
	
done
