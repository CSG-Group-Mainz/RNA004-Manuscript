#this script is a compilation of commands that were run in order to anwser the question how the modification calling for the RNA004 samples had performed, if they had equal throughput to the RNA002 sample
#i e a subsampling of an blood rna004 sample is performed to achieve the same coverage as the rna002 run and then modkit is run as described in the mod_bed_main_samples folder


#load program
ml samtools

#count the read numbers in the RNA002 sample
samtools view -c  RNA002_blood.0.7.2.GRCh38.bam
#5061038
#count the read numbers in the RNA004 sample and express them as a fraction of cooverage
samtools view -c  RNA004_blood.0.7.2.GRCh38.bam
#5061038 /  15017487 = 0,33700964748629381200729522855588
samtools view -c RNA004_S5_DRS_basecall.0.7.2.GRCh38.bam
#5061038 / 17009420 = 0,29754324368496985787875189159889

samtools view -c RNA004_S6_DRS_basecall.0.7.2.GRCh38.bam
#5061038 / 30548388 = 0,1656728335387124191299390331169

#use samtools to subsample to that fraction. 42 serves as seed number
samtools view -@30  -s 42.30 -bh RNA004_S5_DRS_basecall.0.7.2.GRCh38.bam > cov_sub_RNA004_S5_DRS_basecall.0.7.2.GRCh38.bam
samtools view -@30  -s 42.34 -bh RNA004_blood.0.7.2.GRCh38.bam > cov_sub_RNA004_blood.0.7.2.GRCh38.bam
samtools view -@30 -s 42.17 -bh RNA004_S6_DRS_basecall.0.7.2.GRCh38.bam > cov_sub_RNA004_S6_DRS_basecall.0.7.2.GRCh38.bam

#index the subsampled files
samtools index -@8 cov_sub_RNA004_S6_DRS_basecall.0.7.2.GRCh38.bam
samtools index -@8 cov_sub_RNA004_S5_DRS_basecall.0.7.2.GRCh38.bam
samtools index -@8 cov_sub_RNA004_blood.0.7.2.GRCh38.bam

#modkit redo with subsampled bams to have calls with less coverage
#the output of this step will be mod bed files run on the subsampled bams

ml modkit/0.3.1


#step 0 gather samples 

bams_to_use=$( ls cov*38.bam )

echo "start pileup"

function do_pileup {
	  echo $1
		       #
	base=$(basename $1)
	echo $base
	modkit pileup --threads 32 -r GRCh38.primary_assembly.genome.fasta --ignore a  --log-filepath ${base/.bam/.modkit.log}  --filter-threshold T:0.8 --mod-threshold 17802:0.98 $1 ${base/.bam/_pseu.cov_sub_r1.mod.bed}
	modkit pileup --threads 32 --motif DRACH 2 -r GRCh38.primary_assembly.genome.fasta --ignore 17802  --log-filepath ${base/.bam/.modkit.log}  --filter-threshold A:0.8 --mod-threshold a:0.98 $1 ${base/.bam/_m6A.cov_sub_r1.mod.bed}
		}
export -f do_pileup

parallel -v do_pileup ::: $bams_to_use

