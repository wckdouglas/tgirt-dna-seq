#!/usr/bin/env python

PROJECT_PATH=$WORK/cdw2854/ecoli_genome
BED_PATH=$PROJECT_PATH/simulation
BAM_FILE=$PROJECT_PATH/bamFiles
GENOME=$REF/Ecoli/k12_mg1655.fa

#BED_PATH=$WORK/cdw2854/genomeDNA/bamFiles
#BAM_FILE=$BED_PATH
#GENOME=$REF/hg19/Sequence/WholeGenomeFasta/genome.fa

for BED in `ls $BED_PATH/13N*sim*bed | grep -v '[0-9].bed'`
do
	for i in `seq 1 5`
	do
		SAMPLENAME=$(basename ${BED%.bed})
		echo  cat $BED \
			\| bedtools sample -n 2000000 -seed $i -i - \
			\| bedtools nuc -bed - -fi $GENOME -s -seq \
			\| sed 1d \
			\| awk \''{print $4,$NF}'\' OFS=\'\\t\'\
			\| python seq_to_pe_reads.py 75 \
			\| bwa mem -t 4 -p $GENOME - \
			\| samtools sort -@ 4 -O bam -T $BAM_FILE/$SAMPLENAME \
			\> $BAM_FILE/${SAMPLENAME}.${i}.bam
	done
done
