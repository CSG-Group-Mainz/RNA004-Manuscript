#!/usr/bin/env bash

ml samtools
#function get_aln_polys {
# two arguments, a bam file and the tag to extract
#BAM=$1

#TAG=pt
#OUT_PATH=$2

# write a tsv with columns read_id and tag value; we do not want the unaligned reads to show up in the plot lord knows what they are
#echo "$BAM"
#samtools view -@4 -F 4 "$BAM" |  grep -P "pt:i:[0-9]*" |  awk '{print $1"\t"$NF}' | sed "s/pt:i://" > $OUT_PATH/${BAM/.bam/}.poly_a.tsv

#only keep length of unique reads; lets hope so
#sort -s -k1,1 $OUT_PATH/$BAM.poly_a.tsv | uniq > $OUT_PATH/${BAM/.bam/}.poly_a_uniq.tsv

#samtools view -@4 -bh "$BAM" -L /raid/fehof_analysis/RNA004_paper/poly_a_analysis/DDX17/DDX17.bed  > $OUT_PATH/${BAM/.bam/.DDX17.bam}
#samtools view -@4 -bh "$BAM" -L /raid/fehof_analysis/RNA004_paper/poly_a_analysis/OLA1/OLA1.bed  > $OUT_PATH/${BAM/.bam/.OLA1.bam}
#samtools view -@4 -bh "$BAM" -L /raid/fehof_analysis/RNA004_paper/poly_a_analysis/SRP14/SRP14.bed  > $OUT_PATH/${BAM/.bam/.SRP14.bam}
#samtools view -@4 -bh "$BAM" -L /raid/fehof_analysis/RNA004_paper/poly_a_analysis/DDX5/DDX5.bed  > $OUT_PATH/${BAM/.bam/.DDX5.bam}
#samtools view -@4 -bh "$BAM" -L /raid/fehof_analysis/RNA004_paper/poly_a_analysis/RPS24/RPS24.bed  > $OUT_PATH/${BAM/.bam/.RPS24.bam}
#	for i in DDX17 OLA1 SRP14 DX5 RPS24
#		do
#			echo $i
#	samtools view -@4 -F 4 "$OUT_PATH/${BAM/.bam/.$i.bam}" |  grep -P "pt:i:[0-9]*" |  awk '{print $1"\t"$NF}' | sed "s/pt:i://" > $OUT_PATH/${BAM/.bam/$i}.poly_a.tsv
#
#	#only keep length of unique reads; lets hope so
#	sort -s -k1,1 $OUT_PATH/${BAM/.bam/$i}.poly_a.tsv | uniq > $OUT_PATH/${BAM/.bam/$i}.poly_a_uniq.tsv
#		done
#}
#export -f  get_aln_polys
#bams=$( ls ./*38.bam )
#path="/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/DORADO_0_7_2_folder/GLOBAL_polyA_length"
#parallel -v get_aln_polys ::: $bams ::: $path

#function extract_gene_polyA {

# small=$( echo $1 | sed 's/GRCh38//;s/basecall//;s/.bam//' | tr -s "." )
# echo $small
 #we need to use different column matching for different chemistries cause the file sturcture is not the same
# case "$1" in
#	  *"002"*)
 #	samtools view   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$NF}' | sed "s/pt:i:/$small\t/" > ${1/.bam/.poly_a_tags.tsv} ;;
#	  *"004"*)
#	samtools view   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$(NF-3)}' | sed "s/pt:i:/$small\t/" > ${1/.bam/.poly_a_tags.tsv} ;;
#esac
#}
#export -f extract_gene_polyA
#input=$( ls *bam | grep -P -v "HEK|_E"  )
#parallel -v extract_gene_polyA ::: $input
#now get 200000k ivt reads
##function extract200k_IVT {

##	 small=$( basename $1 | sed 's/GRCh38//;s/basecall//;s/.bam//' | tr -s "." | sed -e 's/0.7.2./0.7.2/' )
##	  echo $small
##	case "$1" in
  ##              *"002"*)		
##	  	samtools view   $1 | grep -P "pt:i:[0-9]*" | head -n 400000 | awk '{print $1"\t"$NF}' | sort -u | head -n 200000 | sed "s/pt:i:/$small.0.4.0\t/" > $small.200k.IVT.poly_a_tags.tsv ;;
##		 *"004"*)
##		samtools view   $1 | grep -P "pt:i:[0-9]*" | head -n 400000 | awk '{print $1"\t"$(NF-3)}' | sort -u | head -n 200000 | sed "s/pt:i:/$small.0.4.0\t/" > $small.200k.IVT.poly_a_tags.tsv ;;
##		 esac
 #samtools view   $1 | grep -P "pt:i:[0-9]*" |  awk '{print $1"\t"$NF}' | sort -u  | sed "s/pt:i:/$small\t/" > $small.global.poly_a_tags.tsv
 ##  }
##export -f extract200k_IVT
##input_IVT=$( ls ../*bam | grep "IVT" )
#input_IVT=$( ls ../*bam | grep -P -v "HEK|_E" )
##parallel -v extract200k_IVT ::: $input_IVT

function extract_all {
	small=$( basename $1 | sed 's/GRCh38//;s/basecall//;s/.bam//' | tr -s "." | sed -e 's/0.7.2./0.7.2/' )
	echo $small
	case "$1" in
		*"002"*)
		samtools view -@8 -F 4  $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$NF}' | sed "s/pt:i:/$small\t/" > $small.aligned.global.poly_a_tags.tsv ;;
		*"004"*)
		samtools view -@8 -F 4   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$(NF-3)}' | sort -u | sed "s/pt:i:/$small\t/" > $small.aligned.global.poly_a_tags.tsv ;;
	esac

			 }
export -f extract_all

input_all=$( ls ../*bam | grep -P -v "HEK|_E" )
parallel -v extract_all ::: $input_all



