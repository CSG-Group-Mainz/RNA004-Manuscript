##this script facilitates rebasecalling of all samples with dorado 0.7.2

#load programs (this will only work if your server has this system)
ml samtools/1.16.1
ml minimap2/2.26
ml dorado/0.7.2
ml pod5


#extract sample names and raw data location from file # means do not reprocess samples that have been processed
#the input for this section is a tabular file that lists the sample name in column 1 and the source location of the pod5 data in column 2
for sample in $(cat sample_name_pod5_source.txt | tr "\t" "," | grep -v '#' )
do
bam=$(echo $sample | cut -f1 -d",")   path=$(echo $sample | cut -f2 -d",")
#We need to define different options for either RNA002 samples or RNA004 samples. The automated function in dorado does not work for all of our samples, which is why allocation of options is hardcoded here. But this is a piece of legacy code that should not have to be used for newer dorado versions in combination with newer sequencing runs.
 case "$bam" in
	 *"004"*)
	#echo "$bam   model=RNA004" 
	dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna004_130bps_sup@v5.0.0 $path --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > ${bam/.unaligned.bam/.0.7.2.unaligned.bam} ;;
	*"002"*)
	#echo "$bam   model=RNA002" ;;
	dorado basecaller /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/TEST_folder/rna002_70bps_hac@v3 $path --estimate-poly-a -r --emit-moves --device "cuda:0,cuda:1,cuda:2,cuda:3"  | samtools view -bh - | samtools sort  -@ 32 - -o ${bam/.unaligned.bam/.0.7.2.unaligned.bam} ;;
	*)
	echo "pattern not found" ;;
	 esac
done


#now we align all against grch38
unaligned_bams=$( cat sample_name_pod5_source.txt | tr "\t" "," | grep -v '#' | cut -f1 -d"," | sed 's/.unaligned.bam/.0.7.2.unaligned.bam/' )
echo "align"
##echo $unaligned_bams
##alignment 
function do_align {
samtools fastq -T'*' $1 | minimap2  -y --MD -ax splice -uf -k14 -t 32 /mnt/ssd_share_01/harmonized_rna004_test/references/Homo_sapiens/GENCODE/Sequence/GRCh38.primary_assembly.genome.fa.gz - | samtools sort  -@ 16 - -o ${1/unaligned.bam/GRCh38.bam}
samtools index -@8 ${1/unaligned.bam/GRCh38.bam} 

}
export -f do_align

parallel -v -u --env do_align -j 4 do_align ::: $unaligned_bams

#get bam list and make nicer names for the featureCounts table
bam_list=$( ls *38.bam | tr "\n" " " )
names=$( ls *38.bam | tr "\n" " " | sed 's%/.*/%%' | sed 's/.GRCh38.bam//g; s/_basecall//g' )
#quality control
NanoComp --bam $bam_list -t 40 -n $names --outdir ./reanalysis0.7.2_NanoComp_RNA002_RNA004 --raw --store --tsv_stats
#count features (strandedness has been taken care of via the minimap alignemnt options)
featureCounts -a gencode.v43.basic.annotation.gtf.gz -F 'GTF' -L -s 0 --maxMOp 600 -T 64 -o RNA002_RNA004.DORADO.0.7.2.GRCh38.featurecounts.tsv $bam_list
#source ~/rseqc_env/bin/activate

new=$( ls R*_E*38.bam | tr "\n" "," )
#genebodycoverage
geneBody_coverage.py -r hg38_GENCODE_V45_Basic.bed -i $new  -o E_dorado072.rseqc

