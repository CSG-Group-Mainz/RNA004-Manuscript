#!/usr/bin/env bash
#the purpose of this file is to extract the polyA tags from the aligned bam files and to store them in a tab separated tablr, where column1=unique_read_id and column2=polyAlength (as determined by dorado during the basecalling step)
ml samtools

function extract_all {
	small=$( basename $1 | sed 's/GRCh38//;s/basecall//;s/.bam//' | tr -s "." | sed -e 's/0.7.2./0.7.2/' )
	echo $small
	case "$1" in
		*"002"*)
		samtools view -@8 -F 4  $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$NF}' | sort -u | sed "s/pt:i:/$small\t/" > $small.aligned.global.poly_a_tags.tsv ;;
		*"004"*)
		samtools view -@8 -F 4   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$(NF-3)}' | sort -u | sed "s/pt:i:/$small\t/" > $small.aligned.global.poly_a_tags.tsv ;;
	esac

			 }
export -f extract_all

input_all=$( ls ../*bam | grep -P -v "HEK|_E" )
parallel -v extract_all ::: $input_all


