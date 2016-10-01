#!/bin/bash

DATAPATH=/corral-repl/utexas/2013lambowitz/Data/
TRIMPATH=/corral-repl/utexas/2013lambowitz/Work/Plasma/Trim
PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
TRIMDESTINATION=$PROJECTPATH/trimmedFiles
DESTINATION=$PROJECTPATH/rawData

mkdir -p $DESTINATION

#cat $DATAPATH/JA15513_Yidan/P0311-BZD_S1_R1_001.fastq.gz $datapath/JA15661_Yidan/4-BZD*R1_001.fastq.gz > $DESTINATION/RNA-BZD_R1_001.fq.gz
#cat $DATAPATH/JA15513_Yidan/P0311-BZD_S1_R2_001.fastq.gz $datapath/JA15661_Yidan/4-BZD*R2_001.fastq.gz > $DESTINATION/RNA-BZD_R2_001.fq.gz
#cat $DATAPATH/JA15514_Yidan/P1029-OCD_S2_R1_001.fastq.gz $datapath/JA15661_Yidan/3-OCD-AH_S3_R1_001.fastq.gz > $DESTINATION/RNA-OCD_R1_001.fq.gz
#cat $DATAPATH/JA15513_Yidan/P1029-OCD_S2_R2_001.fastq.gz $datapath/JA15661_Yidan/3-OCD-AH_S3_R2_001.fastq.gz > $DESTINATION/RNA-OCD_R2_001.fq.gz
#cat $DATAPATH/JA15661_Yidan/6-RI*R1_001.fastq.gz > $DESTINATION/RNase-I_R1_001.fq.gz
#cat $DATAPATH/JA15661_Yidan/6-RI*R2_001.fastq.gz > $DESTINATION/RNase-I_R2_001.fq.gz
#cat $TRIMPATH/T2_1.fq.gz $TRIMPATH/T3_1.fq.gz > $DESTINATION/NT_R1_001.fq.gz
#cat $TRIMPATH/T2_2.fq.gz $TRIMPATH/T3_2.fq.gz > $DESTINATION/NT_R2_002.fq.gz

#cp $DATAPATH/JA15930/PD* $DESTINATION

#mv $DESTINATION/SRR2130050_1.fastq.gz $DESTINATION/SRR2130050_R1_001.fq.gz
#mv $DESTINATION/SRR2130050_2.fastq.gz $DESTINATION/SRR2130050_R2_001.fq.gz
#mv $DESTINATION/SRR2130051_1.fastq.gz $DESTINATION/SRR2130051_R1_001.fq.gz
#mv $DESTINATION/SRR2130051_2.fastq.gz $DESTINATION/SRR2130051_R2_001.fq.gz
#mv $DESTINATION/SRR2129993_2.fastq.gz $DESTINATION/SRR2129993_R2_001.fq.gz
#mv $DESTINATION/SRR2129993_1.fastq.gz $DESTINATION/SRR2129993_R1_001.fq.gz

#DATAPATH=/corral-repl/utexas/2013lambowitz/Work/douglas/rawData/jurkatCells/mydata/experiments/preincubation_with_zymol
#cp $DATAPATH/G0-* $DESTINATION
#cp $DATAPATH/G7-* $DESTINATION

# douglas plasma and genome
scp cdw2854@lambcomp01.ccbb.utexas.edu:/stor/work/Lambowitz/Data/NGS/SA16172/dna_plasma_genome/clipped/*gz $DESTINATION

