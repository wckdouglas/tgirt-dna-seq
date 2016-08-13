#!/bin/bash

REFPATH=/scratch/cdw2854/plasmaDNA/reference
GENEFILE=$REFPATH/genes.gtf
OUTPUTFILE=$REFPATH/TSS.gtf

echo Using $GENEFILE
cat $GENEFILE \
    | awk '$3=="gene" && $2=="protein_coding"' \
	> $OUTPUTFILE
echo Written $OUTPUTFILE
