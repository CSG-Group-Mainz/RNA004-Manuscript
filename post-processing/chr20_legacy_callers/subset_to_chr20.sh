#the purpose of this script is to subset the bams to chr 20 and the raw data pod5s. in order to have a small chromosome subset of the data to test various basecalling models on 
#load programs
ml minimap2
ml samtools
#read in data (aligned bams)
bams=$( ls ./DORADO_0_7_2_folder/RNA00*38.bam ) 
#subset bam data to chr20 bam only
function subset {
	
	small_name=$( basename $1 )
	echo $small_name
	samtools view -bh -@16 $1 chr20 > ${small_name/.bam/.chr20.bam}
}
export -f subset
parallel  -u --env subset -j 2 subset ::: $bams

#subset pod5 of A file
mkdir -p chr20_rna002_lemke_a_pod5

#read ids 
#get the read ids out of the bam file to subset the pod5 file to just the subsection of data that aligns to chr20
samtools view RNA002_A_basecall.0.7.2.GRCh38.chr20.bam | cut -f1 | sort -u >read_ids.txt
pod5 filter ./RESEARCH/20230714_Lemke_DRS_RNA002_A_Run2/Sample_A/20230714_1459_2A_PAG68216_8f015b5f/pod5/*.pod5 --output ./chr20_rna002_lemke_a_pod5/chr20_rna002_lemke_a.filtered.pod5 --ids read_ids.txt --missing-ok
mkdir -p chr20_rna002_lemke_a_fast5
#convert to fast5
#because legacy callers take fast5 as input po5 needs to be converted to fast5 format
pod5 convert to_fast5 ./chr20_rna002_lemke_a_pod5/chr20_rna002_lemke_a.filtered.pod5 --output ./chr20_rna002_lemke_a_fast5
#repeat for the other samples
##mkdir -p ./chr20_rna002_lemke_b_pod5/
##samtools view RNA002_B_basecall.0.7.2.GRCh38.chr20.bam | cut -f1 | sort -u >read_idsB.txt
##pod5 filter /mnt/ssd_share_01/new_folder/RESEARCH/20230714_Lemke_DRS_RNA002_B_Run2/Sample_B/20230714_1506_2B_PAG68254_d732e056/pod5/*.pod5 --output ./chr20_rna002_lemke_b_pod5/chr20_rna002_lemke_b.filtered.pod5 --ids read_idsB.txt --missing-ok
##mkdir -p chr20_rna002_lemke_b_fast5
##pod5 convert to_fast5 ./chr20_rna002_lemke_b_pod5/chr20_rna002_lemke_b.filtered.pod5 --output ./chr20_rna002_lemke_b_fast5

##echo "start uhrr"
##mkdir -p ./chr20_rna002_uhrr_pod5/
##samtools view RNA002_UHRR_basecall.0.7.2.GRCh38.chr20.bam | cut -f1 | sort -u >read_idsUHRR.txt
##pod5 filter /mnt/ssd_share_01/new_folder/RESEARCH/20230714_UHRR_R9_DRS/UHRR_R9_DRS/20230714_1508_3A_PAG68321_9f572de7/pod5/*.pod5 --output ./chr20_rna002_uhrr_pod5/chr20_rna002_uhrr_pod5.filtered.pod5 --ids read_idsUHRR.txt --missing-ok
##mkdir -p ./chr20_rna002_uhrr_fast5
##pod5 convert to_fast5 ./chr20_rna002_uhrr_pod5/chr20_rna002_uhrr_pod5.filtered.pod5 --output ./chr20_rna002_uhrr_fast5

