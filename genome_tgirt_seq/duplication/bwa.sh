#!/bin/bash

PROJECTPATH=$WORK/cdw2854/ecoli_genome/UID
DATAPATH=$PROJECTPATH/bam_files/clustered_fastq
RESULPATH=$PROJECTPATH/mapped
INDEX=$REF/Ecoli/b_strain.fa

PROJECTPATH=$WORK/cdw2854/ecoli_genome/rawData/k12
DATAPATH=$PROJECTPATH/umi2id
RESULPATH=$DATAPATH/bam_files
INDEX=$REF/Ecoli/k12_mg1655.fa

THREADS=20
mkdir -p $RESULPATH

for FASTQ1 in $DATAPATH/*R1_001.fastq.gz
do
	SAMPLENAME=$(basename ${FASTQ1%_R1_001.fastq.gz})
	FASTQ2=${FASTQ1/R1/R2}
	echo bwa mem \
		-t $THREADS \
		$INDEX \
		$FASTQ1 \
		$FASTQ2 \
	\| samtools view -@ $THREADS -b \
	\> $RESULPATH/${SAMPLENAME}.bam
done

