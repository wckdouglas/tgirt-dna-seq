#!/bin/bash

WORK_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS

for ANALYTIC in periodicity_gene_body periodicity_tss
do
	for FFT in R scipy
	do
		for FILE in $(ls $WORK_PATH/$ANALYTIC/$FFT/*bed | cut -d'.' -f1 |sort | uniq)
		do
			echo cat $FILE.*.bed \> ${FILE}.bed 
			#echo $FILE
		done
	done
done
