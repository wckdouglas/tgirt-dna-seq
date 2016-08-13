#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/jurkatCells
TABLE_PATH=$PROJECT_PATH/pileup_bases
RESULT_PATH=$PROJECT_PATH/mismatch_table
SNP_BED_FOLDER=/corral-repl/utexas/2013lambowitz/Ref/hg19/annotations/snp_beds
mkdir -p $RESULT_PATH

for TABLE in $TABLE_PATH/*.tsv
do
    SAMPLENAME=$(basename ${TABLE%.tsv})
    CHROM=$(echo ${SAMPLENAME} | rev | cut -d'_' -f1 | rev)
	echo $(which python) mismatch_analysis.py --infile $TABLE \
            --outfile $RESULT_PATH/${SAMPLENAME}.tsv \
            --snp $SNP_BED_FOLDER/${CHROM}.bed
done
