#the purpose of this script is to run mAFiA and m6ABasecaller for blood, UHRR and HEK293T samples
#the input file for mafia is the glori.bed from the nat comm paper for hek293T
#for blood, we used the glori bed that was generated from our in-house sample to serve as input of regions for mAFiA
#it is assumed that the data is already subset to chr20 and converted to fast5, as can be found in subset_to_chr20.sh
ml samtools 
ml guppy/6.2.1

###m6ABasecaller m6A calling
#start basecalling for m6ABasecaller 
guppy_basecaller -i ./chr20_rna002_blood_fast5/ -s ./guppy_m6aBaseCALL_out_blood -c ./m6ABasecaller/basecalling_model/rna_r9.4.1_70bps_m6A_hac.cfg -m ./m6ABasecaller/basecalling_model/template_rna_r9.4.1_70bps_m6A_hac.jsn --fast5_out -x cuda:all

guppy_basecaller -i ./chr20_rna002_blood_IVT_fast5/ -s ./guppy_m6aBaseCALL_out_blood_IVT -c ./m6ABasecaller/basecalling_model/rna_r9.4.1_70bps_m6A_hac.cfg -m ./m6ABasecaller/basecalling_model/template_rna_r9.4.1_70bps_m6A_hac.jsn --fast5_out -x cuda:all

guppy_basecaller -i ./chr20_rna002_lemke_b_fast5/ -s ./guppy_m6aBaseCALL_out_lemke_b -c ./m6ABasecaller/basecalling_model/rna_r9.4.1_70bps_m6A_hac.cfg -m ./m6ABasecaller/basecalling_model/template_rna_r9.4.1_70bps_m6A_hac.jsn --fast5_out -x cuda:all
guppy_basecaller -i ./chr20_rna002_uhrr_fast5/ -s ./guppy_m6aBaseCALL_out_uhrr -c ./m6ABasecaller/basecalling_model/rna_r9.4.1_70bps_m6A_hac.cfg -m ./m6ABasecaller/basecalling_model/template_rna_r9.4.1_70bps_m6A_hac.jsn --fast5_out -x cuda:all


#run m6ABasecaller on data from blood hek293t and uhrr
docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/modphred-3.6.1 /opt/modPhred/run -f /data/GRCh38.primary_assembly.genome.fasta -o /data/m6A_basecaller_blood -i /data/guppy_m6aBaseCALL_out_blood/workspace

docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/modphred-3.6.1 /opt/modPhred/run -f /data/GRCh38.primary_assembly.genome.fasta -o /data/m6A_basecaller_lemke_blood_IVT -i /data/guppy_m6aBaseCALL_out_blood_IVT/workspace

docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/modphred-3.6.1 /opt/modPhred/run -f /data/GRCh38.primary_assembly.genome.fasta -o /data/m6A_basecaller_lemke_b -i /data/guppy_m6aBaseCALL_out_lemke_b/workspace
docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/modphred-3.6.1 /opt/modPhred/run -f /data/GRCh38.primary_assembly.genome.fasta -o /data/m6A_basecaller_lemke_uhrr -i /data/guppy_m6aBaseCALL_out_uhrr/workspace





 ###mAFiA m6A calling
 source ./mAfia/mafia-venv/bin/activate
 

 #download models
 #wget https://zenodo.org/records/8321727/files/models.zip?download=1
 #unpack them
 #unzip 'models.zip?download=1'
 #define varibale paths per sample, here: blood
 models="./RNA004_Manuscript_reanalysis/CHR20_BAMS/RODAN/models"
 output="blood_002_rodan_chr20"
 backbone="${models}/backbone.torch"
 classifiers="${models}/classifiers"
 fast5dir="./chr20_rna002_blood_fast5/"
 ref="GRCh38.primary_assembly.genome.fasta"
 #please note that for blood, the glori bed was genrated from the merged sites of the in house sample as follows:
 #grep "chr20" GLORI_merged_sites.bed | awk -F"\t" '{$6="0"; print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4 }' |  sed  -e '1i chrom\tchromStart\tchromEnd\tname\tscore\tstrand' > glori.bed
 mod="glori.bed"
 mafia="./RNA004_Manuscript_reanalysis/CHR20_BAMS/mAFiA"
 bam="${output}/blood_002_chr20_rodan.q50.bam"

#basecalling with RODAN for mAFiA
#echo "now testing rodan"
python3 ${mafia}/RODAN/basecall.py \
	 --fast5dir ${fast5dir} \
	--model ${backbone} \
	--batchsize 4096 \
	--outdir ${output}

#alignment for mAFiA
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ./RNA004_Manuscript_reanalysis/CHR20_BAMS/UHRR_002_rodan_chr20/rodan.fasta \
	       | samtools view -bST ${ref} -q50 - \
	              | samtools sort - > ${bam}

samtools index ${bam}


#run mAFiA on pre-aligned data
test_mAFiA \
	--bam_file ${bam} \
	--fast5_dir ${fast5dir} \
	--ref_file ${ref} \
	 --mod_file ${mod} \
	--min_coverage 50 \
	--max_num_reads 1000 \
	--backbone_model_path ${backbone} \
	--classifier_model_dir ${classifiers} \
	--mod_prob_thresh 0.5 \
	--out_dir ${output}





