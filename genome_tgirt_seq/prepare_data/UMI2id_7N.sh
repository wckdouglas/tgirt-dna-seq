#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/ecoli_genome/rawData/k12
FASTQ_PATH=${PROJECT_PATH}
OUTPATH=${FASTQ_PATH}/umi2id
PROGRAM_PATH=${HOME}/TGIRT_UMI
PYTHON=$(which python)
mkdir -p ${OUTPATH}

for FQ1 in `ls ${FASTQ_PATH}/K12*7N*R1_001.fastq.gz | grep -v cluster`
do
	FQ2=${FQ1/_R1_/_R2_}
	SAMPLENAME=$(basename ${FQ1%_R1_001.fastq.gz})
	echo ${PYTHON} ${PROGRAM_PATH}/clip_fastq.py \
			--fastq1=${FQ1} \
			--fastq2=${FQ2} \
			--idxBase=13 \
			--constant_region=TAGGATAACAGG \
			--barcodeCutOff=20 \
			--mismatch=7 \
			--outputprefix=- \
		\| python process_7N.py \
		\| ${PYTHON} ${PROGRAM_PATH}/deinterleave_fastq.py \
			- \
			$OUTPATH/${SAMPLENAME}_umi2id_R1_001.fastq.gz \
			$OUTPATH/${SAMPLENAME}_umi2id_R2_001.fastq.gz
done
