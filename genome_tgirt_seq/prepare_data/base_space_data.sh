#!/bin/bash

TARGET_DIR=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12
BASE_SPACE='/stor/scratch/Lambowitz/cdw2854/illumina_basemount/Projects/NextSeq 500 V2: Nextera XT 384-plex E. Coli/Samples/'

for NUM in 5 10 50 60
do
	cat "${BASE_SPACE}"/$NUM/Files/*_R1_001.fastq.gz > $TARGET_DIR/${NUM}_nexteraXT384_R1_001.fastq.gz
	cat "${BASE_SPACE}"/$NUM/Files/*_R2_001.fastq.gz > $TARGET_DIR/${NUM}_nexteraXT384_R2_001.fastq.gz
done
