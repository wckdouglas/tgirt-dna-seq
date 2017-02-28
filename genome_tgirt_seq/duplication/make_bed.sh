#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id
BAM_PATH=$PROJECT_PATH/bam_files
BED_PATH=$BAM_PATH/bed_files
mkdir -p $BED_PATH

for BAM in $BAM_PATH/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo samtools view -bF2048 -F1024 -F512 -F256 -F4 -F8 $BAM \
		\| bedtools bamtobed -mate1 -bedpe \
		\| awk \''$1!="." && $NF!=$(NF-1) && $1==$4 {start=$2;end=$3} {if($5<start) start=$5} {if($6>end) end=$6} {print $1,start,end,$7,end-start,$9}'\' OFS=\'\\t\' \
		\| awk \' '$(NF-1) < 1000' \' \
		\| sed \''s/NC_0/NC-0/g'\' \
		\| tr \'_\' \'\\t\' \
		\| cut -f1,2,3,4,6,7 \
		\| sed \''s/NC-0/NC_0/g'\' \
		\| sort -k1,1 -k2,2n -k3,3n -k6,6 \
		\| datamash -g 1,2,3,6 collapse 4 \
		\| python ~/TGIRT_UMI/dedup_bed/unique_bed.py --infile=- --threshold=3 \
		\> $BED_PATH/${SAMPLENAME}.bed
done
