#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
BED_PATH=$PROJECTPATH/bedFiles/merged/splittedBed
BED_PATH=$PROJECTPATH/bedFiles/merged_rmdup/splittedBed
BED_PATH_UMI=$PROJECTPATH/bedFiles/umi_splitted/demultiplexed
RESULT_PATH=$PROJECTPATH/genomeWPS
REF_PATH=$SCRATCH/GRCh38/hg38_rDNA
GENOME=$REF_PATH/genome_rDNA.fa.fai
PROGRAM=genomeWPS.py
#PROGRAM=strandedGenomeWPS.py
PYTHON=$(which python)

mkdir -p $RESULT_PATH
for BED in  $BED_PATH_UMI/*unique.*.bed #$BED_PATH/*.bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	CHROM=$(echo $SAMPLENAME | rev | cut -d'.' -f1 | rev)
	echo $PYTHON -u $PROGRAM --inFile=$BED \
			--outprefix=$RESULT_PATH/${SAMPLENAME} \
			--genome=$GENOME \
			--TSSwindow=10000 \
			--chromosome=$CHROM
done
