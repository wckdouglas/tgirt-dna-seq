#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/plasmaDNA
FASTQ_PATH=${PROJECT_PATH}/rawData
OUTPATH=${FASTQ_PATH}/umi2id
PROGRAM_PATH=${HOME}/TGIRT_UMI/preproces_fastq
PYTHON=$(which python)
mkdir -p ${OUTPATH}

for FQ1 in `ls ${FASTQ_PATH}/P1*R1_001.fastq.gz | grep -v cluster`
do
	FQ2=${FQ1/_R1_/_R2_}
	SAMPLENAME=$(basename ${FQ1%_R1_001.fastq.gz})
	if [[ $SAMPLENAME == P1022* ]]
	then
		echo ${PYTHON} ${PROGRAM_PATH}/clip_fastq.py \
			--fastq1=${FQ1} \
			--fastq2=${FQ2} \
			--idxBase=13 \
			--barcodeCutOff=20 \
			--constant_region=CATCG \
			--mismatch=2 \
			--outputprefix=- \
		\| ${PYTHON} ${PROGRAM_PATH}/deinterleave_fastq.py \
			- \
			$OUTPATH/${SAMPLENAME}_umi2id_R1_001.fastq.gz \
			$OUTPATH/${SAMPLENAME}_umi2id_R2_001.fastq.gz
	else
		echo ${PYTHON} ${PROGRAM_PATH}/clip_fastq.py \
			--fastq1=${FQ1} \
			--fastq2=${FQ2} \
			--idxBase=13 \
			--barcodeCutOff=20 \
			--outputprefix=- \
		\| ${PYTHON} ${PROGRAM_PATH}/deinterleave_fastq.py \
			- \
			$OUTPATH/${SAMPLENAME}_umi2id_R1_001.fastq.gz \
			$OUTPATH/${SAMPLENAME}_umi2id_R2_001.fastq.gz
	fi
done

