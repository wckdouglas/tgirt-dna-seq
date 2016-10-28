for PROFILE in profiles/P1203-SQ1_S2.csv profiles/SRR2130051.csv
do
	if [[ $PROFILE == *SRR* ]]
	then
		SAMPLENAME=ssDNA-seq_sim
	else
		SAMPLENAME=TGIRT_sim
	fi
	echo python tgirt_seq_simulator.py \
		-i profiles/insert_profile.csv \
		-b $PROFILE \
		-r $REF/GRCh38/hg38_rDNA/genome_rDNA.fa \
		-f 10 \
		-o  $WORK/cdw2854/plasmaDNA/bedFiles/$SAMPLENAME \
		-t 23
done
