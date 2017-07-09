#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/bisufite_seq
RAW_DATA=$PROJECT_PATH/raw_Data
OUT_PATH=$RAW_DATA/UMI_clustered
PROGRAM=/stor/home/cdw2854/TGIRT_UMI/genome_map_reduce/meth_R1R_13N_pipeline.py 
INDEX=$REF/hg19/Sequence/WholeGenomeFasta/genome.fa
THREADS=8

for FQ1 in $(ls $RAW_DATA/*R1_001.fastq.gz | grep 'PB\|P13B')
do
	FQ2=${FQ1/R1_001/R2_001}
	SAMPLENAME=$(basename ${FQ1%_R1_001.fastq.gz})
	echo $(which python) ${PROGRAM} \
		-1 ${FQ1} -2 ${FQ2} -o $OUT_PATH \
		-x $INDEX -t $THREADS \
		\&\> $OUT_PATH/${SAMPLENAME}.log
done
