#!/bin/bash

python extractFpkmTss.py \
    /scratch/02727/cdw2854/plasma_project/TGIRT_plasma_dna/assemblies/RNA-OCD.gtf \
    | uniq \
    | sort -k8nr \
    | awk '$8>0 && $2>0' \
    > /scratch/02727/cdw2854/plasma_project/TGIRT_plasma_dna/assemblies/genes.tss.bed
