#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/target-seq/raw_reads
DATA_PATH=$PROJECT_PATH
RESULT_PATH=$PROJECT_PATH/splitted
SUFFIX=_R1_001.fastq.gz
PROGRAM_PATH=$HOME/indexed-tgirt-seq
PROGRAM=$PROGRAM_PATH/read_cluster_pairs.py
THREADS=12
mkdir -p  $RESULT_PATH

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
	echo $(which python)  $PROGRAM \
		--outputprefix ${RESULT_PATH}/${SAMPLE_NAME}-errorFree \
	    --fastq1 ${FQ1} \
		--fastq2 ${FQ2} \
		--idxBase 15 \
		--barcodeCutOff 20 \
	    --cutoff 0 \
		--constant_region TACGCTCTTTCTCCGCGAATGCGGCGAGCGAGCTGGATGTC \
		--threads $THREADS \
		--mismatch 3 \
		--read read2 \
		--fraction 0.66 \
		\&\> ${RESULT_PATH}/${SAMPLE_NAME}.log
done

#PROGRAM=double_index_cluster.py
#for FQ1 in `ls $DATA_PATH/DB*${SUFFIX}`
#do
#	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
#	FQ2=${FQ1/R1/R2}
#	echo $(which python) $PROGRAM \
#		--outputprefix=$RESULT_PATH/$SAMPLE_NAME-errorFree-double-BC \
#	    --fastq1=$FQ1 --fastq2=$FQ2 \
#		--idxBase=13 --barcodeCutOff=30 \
#	    --cutoff 0 -l CATCG -r GAGTGTAGTGCATATGAGCACTGTCGAT \
#		--threads $THREADS
#done
