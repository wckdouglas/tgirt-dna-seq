#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/ecoli_genome/
BAM_PATH=$PROJECT_PATH/picard_results
OUT_PATH=$PROJECT_PATH/fragment_ends
PROGRAM=/stor/home/cdw2854/tgirt-dna-seq/genome_tgirt_seq/simulator/profiles/seq_to_basecount.py
for BAM in $BAM_PATH/*.MarkDuplicate.bam 
do 
	SAMPLENAME=$(basename ${BAM%.bam})
	echo samtools sort -n -T $OUT_PATH/$SAMPLENAME -O bam $BAM \
		\| bedtools bamtobed -i - -bedpe -mate1  \
		\| bedpe_to_bed.py \
		\| bedtools nuc -bed - -fi $REF/Ecoli/k12_mg1655.fa -s -seq \
		\| sed 1d \| rev \| cut -f1 \| rev \
		\| python $PROGRAM $OUT_PATH/$SAMPLENAME
done

