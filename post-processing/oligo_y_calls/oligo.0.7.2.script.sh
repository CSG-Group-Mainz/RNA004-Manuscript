##this script facilitates rebasecalling of all lemke oligos rna004 samples with dorado 0.7.2
##input data is ca 100.00s reads from felix folder
ml samtools/1.16.1
ml minimap2/2.26
ml dorado/0.7.2
ml pod5

#now we basecall the pod5 from felix folder

dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/fehof_analysis/RNA004_paper/nanoCEM/Oligo1/ --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > Oligo1.0.7.2.unaligned.bam 

dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/fehof_analysis/RNA004_paper/nanoCEM/Oligo2/ --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > Oligo2.0.7.2.unaligned.bam


#now we check what happens if we mess with the modified bases model and force it to evalute Y whether that fixes C/T calls or not
#now we basecall the pod5 from felix folder
dorado download --model rna004_130bps_hac@v5.0.0_m6A@v1
#cd there and change config.toml to Y
dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/fehof_analysis/RNA004_paper/nanoCEM/Oligo1/  -r --emit-moves --modified-bases-models ./rna004_130bps_sup@v5.0.0_pseU@v1 -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > Y_eval.Oligo1.0.7.2.unaligned.bam

dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/fehof_analysis/RNA004_paper/nanoCEM/Oligo2/  -r --emit-moves --modified-bases-models ./rna004_130bps_sup@v5.0.0_pseU@v1 -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > Y_eval.Oligo2.0.7.2.unaligned.bam


#and now once again for the vektor samples
dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/awiercze_analysis/RESEARCH/LEMKE_SAMPLES/RNA004_A.ROIv5.pod5  --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > A_vector_004.0.7.2.unaligned.bam

dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/awiercze_analysis/RESEARCH/LEMKE_SAMPLES/RNA004_B.ROIv5.pod5  --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > B_vector_004.0.7.2.unaligned.bam


#now we check what happens if we  force the pseu models to evaluate miscalled c bases. whether modprobs are okay or not
#now we basecall the pod5 from felix folder
##dorado download --model rna004_130bps_hac@v5.0.0_m6A@v1
#cd there and change config.toml to Y
dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/awiercze_analysis/RESEARCH/LEMKE_SAMPLES/RNA004_A.ROIv5.pod5   -r --emit-moves --modified-bases-models ./rna004_130bps_sup@v5.0.0_pseU@v1 -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > Y_eval.A_vector_004.0.7.2.unaligned.bam

dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 /raid/awiercze_analysis/RESEARCH/LEMKE_SAMPLES/RNA004_B.ROIv5.pod5   -r --emit-moves --modified-bases-models ./rna004_130bps_sup@v5.0.0_pseU@v1 -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > Y_eval.B_vector_004.0.7.2.unaligned.bam
#########



#now we align all against oligo/vector reference
unaligned_bams=$( ls *vector*.0.7.2.unaligned.bam )
echo "align"
echo $unaligned_bams
#alignment 
function do_align {
samtools fastq -T'*' $1 | minimap2  -y --MD -ax splice -uf -k14 -t 32 reference_LUKAS.fa - | samtools sort  -@ 16 - -o ${1/unaligned.bam/vector.bam}
samtools index -@8 ${1/unaligned.bam/vector.bam} 

}
export -f do_align

parallel -v -u --env do_align -j 4 do_align ::: $unaligned_bams

#get modbeds from modbams
ml modkit

#do get basemods but just for rna004 since its not meaningful for rna002
bams=$( ls *oligo.bam )
function do_get_modbeds {
	        small_name=$( basename $1 )
		modkit pileup -t 4 $1 ${small_name/.bam/.mod.bed}

		}
export -f do_get_modbeds
parallel -v -u --env do_get_modbeds -j 4 do_get_modbeds ::: $bams


#modkit sample probs to get the histogram of modification probabilities
modkit sample-probs -t 32 --log-filepath oligo1.log  -p 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 --prefix oligo1 --hist --only-mapped Oligo1.0.7.2.oligo.bam --out-dir sample_probs_oligo1

