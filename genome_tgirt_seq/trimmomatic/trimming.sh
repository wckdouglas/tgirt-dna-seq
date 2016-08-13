#!/bin/bash


PROJECT_PATH=/scratch/02727/cdw2854/jurkatCells
FASTQ_PATH=$PROJECT_PATH/splitted
TRIM_PATH=$PROJECT_PATH/trimmed
THREADS=6
SUFFIX=_R1_001.fastq.gz
mkdir -p $TRIM_PATH

for FQ in $FASTQ_PATH/*$SUFFIX
do
	SAMPLENAME=$(basename ${FQ%$SUFFIX})
	if [[ $SAMPLENAME == *DB* ]]
	then
		ADAPTOR=double_index_adaptors.fa
	elif [[ $SAMPLENAME == *SB* ]]
	then
		ADAPTOR=indexed_adaptors.fa
	else
		echo specify adaptor?
		exit
	fi
	OPTIONS="ILLUMINACLIP:$ADAPTOR:1:20:20:6:true MINLEN:10"
	echo trimmomatic PE -threads $THREADS \
		-basein $FQ -baseout $TRIM_PATH/${SAMPLENAME}.fq.gz \
		$OPTIONS
done
