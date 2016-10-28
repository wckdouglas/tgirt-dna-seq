echo python tgirt_seq_simulator.py \
	-i profiles/insert_profile.csv \
	-b profiles/dinucleotide_profile.csv \
	-r $REF/GRCh38/hg38_rDNA/genome_rDNA.fa \
	-f 10 \
	-o  $WORK/cdw2854/plasmaDNA/bedFiles/tgirt_sim
