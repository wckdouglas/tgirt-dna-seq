for SIDE in 3 5 both no
do
	FOLD=1000
	for PROFILE_PREFIX in 13N_clustered
	do
		if [[ $SIDE == both ]]
		then
			SAMPLENAME=K12_sim
			FOLD=10000
		elif [[ $SIDE == 5 ]]
		then
			SAMPLENAME=K12_sim_ligation_only
		elif [[ $SIDE == no ]]
		then
			SAMPLENAME=K12_sim_no_bias
			FOLD=10
		elif [[ $SIDE == 3 ]]
		then
			SAMPLENAME=K12_sim_template_switch
		fi
		
		# E coli
		#	-i ./profiles/13N_len_profile.csv \
		#	-b ./profiles/13N_base_profile.csv \
		echo python tgirt_seq_kmer_simulator.py \
			-r $REF/Ecoli/k12_mg1655.fa \
			-f $FOLD \
			-i ./profiles/${PROFILE_PREFIX}_len_profile.csv \
			-b ./profiles/${PROFILE_PREFIX}_base_profile.csv \
			-o $WORK/cdw2854/ecoli_genome/simulation/${PROFILE_PREFIX}_${SAMPLENAME} \
			-t 20 \
			-s $SIDE \
			\&\> log/${PROFILE_PREFIX}_${SAMPLENAME}.log
	done
done
