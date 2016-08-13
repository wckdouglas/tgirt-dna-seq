#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
PROJECTPATH=/stor/work/Lambowitz/cdw2854/plasmaDNA
BED_PATH=$PROJECTPATH/rmdup_bed/splittedBed
RESULT_PATH=$PROJECTPATH/genomeWPS
REF_PATH=/corral-repl/utexas/2013lambowitz/Ref/GRCh38/hg38_rDNA
REF_PATH=/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA
GENOME=$REF_PATH/genome_rDNA.fa.fai
PROGRAM=genomeWPS.py
#PROGRAM=strandedGenomeWPS.py
PYTHON=$(which python)

mkdir -p $RESULT_PATH
for BED in $BED_PATH/*bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	CHROM=$(echo $SAMPLENAME | rev | cut -d'.' -f1 | rev)
	echo $PYTHON -u $PROGRAM --inFile=$BED \
			--outprefix=$RESULT_PATH/${SAMPLENAME} \
			--genome=$GENOME \
			--TSSwindow=10000 \
			--chromosome=$CHROM
done
