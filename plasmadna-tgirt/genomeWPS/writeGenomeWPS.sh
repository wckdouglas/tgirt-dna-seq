#!/bin/bash

PROJECTPATH=$WORK/cdw2854/plasmaDNA
PROJECTPATH=$SCRATCH/plasmaDNA
BED_PATH=$PROJECTPATH/bedFiles/splittedBed
BED_PATH=$PROJECTPATH/bedFiles/chrom_split_bed/demultiplexed
RESULT_PATH=$PROJECTPATH/genomeWPS
REF_PATH=$REF/GRCh38/hg38_rDNA
GENOME=$REF_PATH/genome_rDNA.fa.fai
PROGRAM=genomeWPS.py
#PROGRAM=strandedGenomeWPS.py
PYTHON=$(which python)

mkdir -p $RESULT_PATH
for BED in  $BED_PATH/*unique.*.bed #$BED_PATH/*.bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	CHROM=$(echo $SAMPLENAME | rev | cut -d'.' -f1 | rev)
	echo $PYTHON -u $PROGRAM --inFile=$BED \
			--outprefix=$RESULT_PATH/${SAMPLENAME} \
			--genome=$GENOME \
			--TSSwindow=10000 \
			--chromosome=$CHROM
done
