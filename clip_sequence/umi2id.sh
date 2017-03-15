#!/bin/bash

DATAPATH=$SCRATCH/plasmaDNA/rawData
CLIPPED_PATH=$DATAPATH/umi2id
SOFTWARE_PATH=$HOME/TGIRT_UMI/preproces_fastq 
mkdir -p $CLIPPED_PATH

for FQ1 in $DATAPATH/*R1_001.fastq.gz
do
	R1_NAME=$(basename $FQ1)
	FQ2=${FQ1/R1/R2}
	R2_NAME=$(basename $FQ2)
	echo python $SOFTWARE_PATH/clip_fastq.py \
			--fastq1 $FQ1 \
			--fastq2 $FQ2 \
			--idxBase 13 \
			--barcodeCutOff 20  \
		\| python $SOFTWARE_PATH/deinterleave_fastq.py - \
			$CLIPPED_PATH/${R1_NAME%_R1_001.fastq.gz}_umi2id_R1_001.fastq.gz \
			$CLIPPED_PATH/${R1_NAME%_R2_001.fastq.gz}_umi2id_R2_001.fastq.gz 
done
