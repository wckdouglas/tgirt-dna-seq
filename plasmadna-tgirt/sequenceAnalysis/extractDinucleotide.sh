#!/bin/bash

PROJECTPATH=/Users/wckdouglas/plasmaDNA/results
BEDPATH=$PROJECTPATH/bedFile
RESULTPATH=$PROJECTPATH/dinucleotides
REFFASTA=/Users/wckdouglas/plasmaDNA/reference/genome_rDNA.fa
PYTHON=$(which python)
PROGRAM=dinucleotide.py
mkdir -p $RESULTPATH

for BED in $BEDPATH/*bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	echo $PYTHON $PROGRAM --inBed=$BED \
		--outprefix=$RESULTPATH/$SAMPLENAME \
		--refFasta=$REFFASTA \
		--window=400
done
