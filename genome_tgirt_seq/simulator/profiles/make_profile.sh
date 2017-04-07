BAM_FILE=/stor/work/Lambowitz/cdw2854/ecoli_genome/picard_results/K12_UMI_1_S9.MarkDuplicate.bam 
BAM_FILE=/stor/work/Lambowitz/cdw2854/ecoli_genome/picard_results/K12_F_NEB_S7_umi2id.MarkDuplicate.bam
BAM_FILE=/stor/work/Lambowitz/cdw2854/ecoli_genome/bamFiles/mark_duplicate/K12_UMI_1_S9_clustered.bam

samtools sort -n $BAM_FILE \
	| bamToBed -bedpe -mate1 \
	| awk '$1 != "\."' \
	| bedpe_to_bed.py \
	| bedtools nuc -s -fi $REF/Ecoli/k12_mg1655.fa -bed - -seq \
	| awk '{print $NF}' \
	| python seq_to_kmer.py 13N_clustered 8

