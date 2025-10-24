#### CHECK JOINED COVERED SITES 
#!/bin/bash
ml samtools
# Requirements: samtools, GNU parallel
# Input: bam_list.txt and bed_list.txt must exist

### SUBSAMPLED BLOOD
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/R2Q2/cov_sub_RNA004_*.bam > /home/awiercze/BAM_FILE_LIST_subsampled_blood.txt
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/DORADO_0_7_2_folder/RNA002_blood.0.7.2.GRCh38.bam >> /home/awiercze/BAM_FILE_LIST_subsampled_blood.txt
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/CHR20_BAMS/RNA002_folder/m6A_*blood*.bed > /home/awiercze/BED_FILE_LIST_subsampled_blood.txt
ls /home/awiercze/final_blood_002_rodan_chr20_mAFiA.sites.headerremoved.bed >> /home/awiercze/BED_FILE_LIST_subsampled_blood.txt
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/R2Q2/cov_sub_RNA004*bed >> /home/awiercze/BED_FILE_LIST_subsampled_blood.txt

### ALL M6A
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/m6A_pseu_refilter_rev1_nar1/*m6A*.bed > /home/awiercze/BED_FILE_LIST.txt
ls /home/awiercze/final_blood_*_002_rodan_chr20_mAFiA.sites.headerremoved.bed >> /home/awiercze/BED_FILE_LIST.txt
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/CHR20_BAMS/RNA002_folder/m6A_*blood*.bed >> /home/awiercze/BED_FILE_LIST.txt

### ALL PSEU
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/m6A_pseu_refilter_rev1_nar1/*pseu*.bed > /home/awiercze/PSEU_BED_FILE_LIST.txt
ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S6_IVT_analysis/merged_RNA004_*_pseu.r1.mod.bed >> /home/awiercze/PSEU_BED_FILE_LIST.txt


# Create a combined list with 3 columns: BAM, BED, 
: '
for i in $(cat /home/awiercze/PSEU_BED_FILE_LIST.txt); do
    echo "$i"
    baseName=$(basename $i)
    samtools depth -aa -b $i -@ 100 -f BAM_FILE_LIST_004.txt -H -o /home/awiercze/${baseName/.bed/.SAMPLES_DEPTH.txt}
done
'

: '
for i in $(cat /home/awiercze/BED_FILE_LIST.txt); do
    echo "$i"
    baseName=$(basename $i)
    samtools depth -aa -b $i -@ 100 -f BAM_FILE_LIST_004.txt -H -o /home/awiercze/${baseName/.bed/.SAMPLES_DEPTH.txt}
done
'

for j in $(cat /home/awiercze/BED_FILE_LIST_subsampled_blood.txt); do
    for k in $(cat /home/awiercze/BAM_FILE_LIST_subsampled_blood.txt); do
        bam_name=$( basename $k)
        bam_name=${bam_name/.bam/}
        echo "$j"
        echo $k
        wc -l $j
        baseName=$( basename $j )
        samtools depth -aa -b $j -@ 100 -H -o /home/awiercze/${baseName/.bed/.SUBSAMPLED_DEPTH}_$bam_name.txt $k
    done
   
done