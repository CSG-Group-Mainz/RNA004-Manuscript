#this script is a compilation of commands that are needed to reprocess the ont epigenetic data
#download pod5 references and bam data from the aws repo
aws s3 sync --no-sign-request s3://ont-open-data/rna-modbase-validation_2025.03/references references
aws s3 sync --no-sign-request s3://ont-open-data/rna-modbase-validation_2025.03/subset subset
aws s3 sync --no-sign-request s3://ont-open-data/rna-modbase-validation_2025.03/basecalls test
#re-basecall pod5 data and align
ml samtools/1.16.1
ml minimap2/2.26
ml dorado/0.7.2
ml modkit/0.3.1
##echo "step1"
dorado basecaller ./rna004_130bps_sup@v5.0.0 ./ONT_EPI_RNA_TEST_DATA/subset/m6A_rep1.pod5  --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > RNA004_ONT.M6A.R1.unaligned.bam
dorado basecaller ./rna004_130bps_sup@v5.0.0 ./ONT_EPI_RNA_TEST_DATA/subset/m6A_rep2.pod5  --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > RNA004_ONT.M6A.R2.unaligned.bam


dorado basecaller ./rna004_130bps_sup@v5.0.0 ./pseU_rep1.pod5  --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > RNA004_ONT.pseu.R1.unaligned.bam
dorado basecaller ./rna004_130bps_sup@v5.0.0 ./pseU_rep2.pod5  --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > RNA004_ONT.pseu.R2.unaligned.bam


echo "step 2 align"
samtools fastq -T'*' RNA004_ONT.M6A.R1.unaligned.bam | minimap2  -y --MD -ax splice -uf -k14 -t 32 ./references/sampled_context_strands.fa - | samtools sort  -@ 16 - -o RNA004_ONT.M6A.R1.EPIref.bam
samtools fastq -T'*' RNA004_ONT.M6A.R2.unaligned.bam | minimap2  -y --MD -ax splice -uf -k14 -t 32 ./references/sampled_context_strands.fa - | samtools sort  -@ 16 - -o RNA004_ONT.M6A.R2.EPIref.bam


samtools fastq -T'*' RNA004_ONT.pseu.R1.unaligned.bam | minimap2  -y --MD -ax splice -uf -k14 -t 32 ./references/sampled_context_strands.fa - | samtools sort  -@ 16 - -o RNA004_ONT.pseu.R1.EPIref.bam
samtools fastq -T'*' RNA004_ONT.pseu.R2.unaligned.bam | minimap2  -y --MD -ax splice -uf -k14 -t 32 ./references/sampled_context_strands.fa - | samtools sort  -@ 16 - -o RNA004_ONT.pseu.R2.EPIref.bam

echo "step3 modkit"
bams_to_use=$( ls  *ref.bam)

function do_pileup {
	 echo $1
	samtools index -@16 $1
	base=$1
	modkit sample-probs --threads 8 $1 --hist -o hist_dir${base/.bam/}  --percentiles 0.1,0.8,0.85,0.9
}
export -f do_pileup
parallel -v do_pileup ::: $bams_to_use
