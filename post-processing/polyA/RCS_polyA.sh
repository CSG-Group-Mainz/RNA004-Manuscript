ml samtools/1.16.1
ml minimap2/2.26
ml dorado/0.7.2
#this script is written to compare poly(A) lengths of a standardized sequence (here RCS -the reference control strand that is in every ONT kit) between RNA004 and RNA002. The RCS was sequenced as a spike-in for the first runs of uhrr

unaligned_bams=$( ls *UHRR*unaligned.bam )

echo "align"
echo $unaligned_bams
#alignment
function do_align {
#only output aligned hits cause we dont have the space (makes bams smaller)
samtools fastq -T'*' $1 | minimap2  -y --MD -ax splice -uf --sam-hit-only -k14 -t 32 Enolase_RCS_ONT.fasta - | samtools sort  -@ 16 - -o ${1/unaligned.bam/RCS.aligned.bam}

samtools index -@8 ${1/unaligned.bam/.RCS.aligned.bam}

}

#export -f do_align
parallel -v do_align ::: $unaligned_bams


#get the polyAtags in a text file
echo "now extract"
function extract {

      
       small=$( basename $1 | sed 's/.bam//'  )
        echo $small
      case "$1" in
                *"002"*)
                samtools view   $1 | grep -P "pt:i:[0-9]*"  | awk '{print $1"\t"$NF}' | sort -u | sed "s/pt:i:/$small\t/" > $small.poly_a_tags.tsv ;;
               *"004"*)
	              samtools view   $1 | grep -P "pt:i:[0-9]*" | awk '{print $1"\t"$(NF-3)}' | sort -u | sed "s/pt:i:/$small\t/" > $small.poly_a_tags.tsv ;;
      esac             
           }
 export -f extract

input_IVT=$( ls *RCS.bam )
 
 parallel -v extract ::: $input_IVT

 
