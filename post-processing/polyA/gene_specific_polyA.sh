#!/usr/bin/env bash
#this script extract gene specific poly(A) lengths
#you will need a BED file per gene and an aligned bam as input
#output is a table: with read_id <TAB> poly(Alength) of read
ml samtools
function get_aln_polys {
# two arguments, a bam file and the tag to extract
BAM=$1

TAG=pt
OUT_PATH=$2


#subset the sample bam to smaller bam files per gene
samtools view -@4 -bh "$BAM" -L ./poly_a_analysis/DDX17/DDX17.bed  > $OUT_PATH/${BAM/.bam/.DDX17.bam}
samtools view -@4 -bh "$BAM" -L ./poly_a_analysis/OLA1/OLA1.bed  > $OUT_PATH/${BAM/.bam/.OLA1.bam}
samtools view -@4 -bh "$BAM" -L ./poly_a_analysis/SRP14/SRP14.bed  > $OUT_PATH/${BAM/.bam/.SRP14.bam}
samtools view -@4 -bh "$BAM" -L ./poly_a_analysis/DDX5/DDX5.bed  > $OUT_PATH/${BAM/.bam/.DDX5.bam}
samtools view -@4 -bh "$BAM" -L ./poly_a_analysis/RPS24/RPS24.bed  > $OUT_PATH/${BAM/.bam/.RPS24.bam}

}
export -f  get_aln_polys
bams=$( ls ./*38.bam )
path="/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/DORADO_0_7_2_folder/GLOBAL_polyA_length"
parallel -v get_aln_polys ::: $bams ::: $path

#now the polya tags are extracted from the subset bams and stored in a table
function extract_gene_polyA {

 small=$( echo $1 | sed 's/GRCh38//;s/basecall//;s/.bam//' | tr -s "." )
 echo $small
 #we need to use different column matching for different chemistries cause the file sturcture is not the same
 case "$1" in
         *"002"*)
       samtools view   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$NF}' | sed "s/pt:i:/$small\t/" > ${1/.bam/.poly_a_tags.tsv} ;;
         *"004"*)
       samtools view   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$(NF-3)}' | sed "s/pt:i:/$small\t/" > ${1/.bam/.poly_a_tags.tsv} ;;
esac
}
export -f extract_gene_polyA

input=$( ls *bam )
parallel -v extract_gene_polyA ::: $input




