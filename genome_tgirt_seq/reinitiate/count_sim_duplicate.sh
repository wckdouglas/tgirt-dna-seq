#!/bin/bash

# K12_UMI_1_S9_umi2id.unique.bed
cat /stor/work/Lambowitz/cdw2854/ecoli_genome/simulation/13N_clustered_K12_sim.bed \
	| head -n 2098603 \
	| cut -f1-3 \
	| uniq -c  \
	| awk '{print $1}'  \
	| sort \
	| uniq -c 

