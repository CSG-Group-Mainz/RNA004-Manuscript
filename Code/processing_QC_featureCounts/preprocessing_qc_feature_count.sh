#_____________________________________________________________________________________________
# RNA002 vs RNA004 SCRIPT: Feature Quantification and Quality Checks


# Load tools
ml pod5
ml dorado
ml samtools
ml minimap2
ml bcftools

# Basecalling for RNA002
function run_dorad_RNA002 {
	time dorado basecaller \
        ./rna002_70bps_hac@v3 $1 \
        --estimate-poly-a \
        -r \
        --emit-moves \
        --device "cuda:0,cuda:1,cuda:2,cuda:3" \
        | samtools view -bh - \
        | samtools sort  -@ 32 - -o "$2"_basecall.unaligned.bam
}
export -f run_dorad_RNA002

# Basecalling for RNA004
function run_dorado_RNA004 {
	time dorado basecaller \
        ./rna004_130bps_sup@v3.0.1 $1 \
        --modified-bases m6A_DRACH \
        --estimate-poly-a \
        -r \
        --device "cuda:0,cuda:1,cuda:2,cuda:3" \
        | samtools view -bh - \
        | samtools sort  -@ 32 - -o "$2"_basecall.unaligned.bam
}
export -f run_dorado_RNA004

# Mapping onto GChr38
function mapping_GChr38 {
	samtools fastq -T '*' "$1" \
	| minimap2 -y --MD -ax splice -uf -k14 -t 96 $2 - \
	| samtools view -bh \
	| samtools sort -@32 - -o ${1/.bam/}.GChr38.bam
	samtools index ${1/.bam/}.$2.bam
}
export mapping_GRChr38

# RSeQCs genebodycoverage and TIN function
function run_rseqc {
    geneBody_coverage.py -r /home/awiercze/hg38_GENCODE_V45_Basic.bed -i $1  -o ${1/.bam}
    tin.py -r /home/awiercze/hg38_GENCODE_V45_Basic.bed -i $1 
}
export -f run_rseqc

# GeneBodyCoverage for single gene 
function get_single_gene_cov {
    gene=$2
    positions=$(zgrep $gene gencode.v43.annotation.gtf.gz | awk '$3=="gene"' | cut -f1,4,5 | sed -e 's/\t/:/;s/\t/-/')
    echo $1
    samtools view -bh $1 "$positions" > ${1/.bam/}.$gene.bam
    samtools index ${1/.bam/}.$gene.bam
    geneBody_coverage.py -r hg38_GENCODE_V45_Basic.bed -i ${1/.bam/}.$gene.bam -o ${1/.bam/}.$gene
}
export -f get_single_gene_cov


#____________________________________________________________________________________________________
# BASECALLING, MAPPING AND QUALITY CONTROL OF ALL SAMPLES

# Run Dorado for blood samples 
RNA004_A="./RNA004_HEK293_A_pod5/"
RNA004_B="./RNA004_HEK293_B_pod5/"
RNA002_A="./RNA002_HEK293_A_pod5/"
RNA002_B="./RNA002_HEK293_B_pod5/"

RNA002_UHRR="./RNA002_UHRR_pod5/"
RNA004_UHRR_1="./RNA004_UHRR_1_pod5/"
RNA004_UHRR_2="./RNA004_UHRR_2_pod5/"

RNA004_IVT="RNA004_runs/RNA004_peripheral_blood/IVT"
RNA004_Direct="RNA004_runs/RNA004_peripheral_blood/Normal"

run_dorad_RNA002 $RNA002_A "RNA002_A"
run_dorad_RNA002 $RNA002_B "RNA002_B"
run_dorado_RNA004 $RNA004_A "RNA004_A"
run_dorado_RNA004 $RNA004_B "RNA004_B"
run_dorad_RNA002 $RNA002_UHRR "RNA002_UHRR"
run_dorado_RNA004 $RNA004_UHRR_1 "RNA004_UHRR_1"
run_dorado_RNA004 $RNA004_UHRR_2 "RNA004_UHRR_2"
run_dorado_RNA004 $RNA004_IVT "RNA004_blood_IVT"
run_dorado_RNA004 $RNA004_Direct "RNA004_blood_Direct"

# Map all Samples to complete GRCh38 genome (GENCODE Release 43)
RNA002_RNA004_unaligned_bams=$( ls *unaligned.bam )
parallel -v -u --env mapping_GChr38 -j 9 mapping_GChr38 ::: $RNA002_RNA004_unaligned_bams ::: GRCh38.primary_assembly.genome.fa


# Extract global genebodycoverage and TIN values (using GENCODE bed file from RSEQC)
RNA002_RNA004_aligned_bams=$( ls *GRCh38.bam )

parallel -v -u --env run_rseqc -j 9 run_rseqc ::: $RNA002_RNA004_aligned_bams

# Run genebodycoverage on single gene level for all samples (using GENCODE gtf Release 43) 

parallel -v -u --env get_single_gene_cov -j 9 get_single_gene_cov ::: $RNA002_RNA004_aligned_bams ::: XIST 
parallel -v -u --env get_single_gene_cov -j 9 get_single_gene_cov ::: $RNA002_RNA004_aligned_bams ::: MT-ND3 

# Check Quality using NanoComp
NanoComp --bam \
    RNA002_A_basecall.GRCh38.bam RNA002_B_basecall.GRCh38.bam \
    RNA002_UHRR_basecall.GRCh38.bam \
    RNA004_A_basecall.GRCh38.bam RNA004_B_basecall.GRCh38.bam \
    RNA004_UHRR_1_basecall.GRCh38.bam RNA004_UHRR_2_basecall.GRCh38.bam \
    RNA004_blood_Direct_basecall.GRCh38.bam RNA004_blood_IVT_basecall.GRCh38.bam \
    -t 40 \
    -n RNA002_HEK293_A RNA002_HEK293_B RNA002_UHRR RNA004_HEK293_A RNA004_HEK293_B RNA004_UHRR_1 RNA004_UHRR_2 RNA004_blood_normal RNA004_blood_IVT \
    --outdir /home/awiercze/NanoComp_RNA002_RNA004 --raw --store --tsv_stats


NanoComp --ubam \
    RNA002_A_basecall.unaligned.bam RNA002_B_basecall.unaligned.bam \
    RNA002_UHRR_basecall.unaligned.bam \
    RNA004_A_basecall.unaligned.bam RNA004_B_basecall.unaligned.bam \
    RNA004_UHRR_1_basecall.unaligned.bam RNA004_UHRR_2_basecall.unaligned.bam \
    RNA004_blood_Direct_basecall.unaligned.bam RNA004_blood_IVT_basecall.unaligned.bam \
    -t 94 \
    -n RNA002_HEK293_A RNA002_HEK293_B RNA002_UHRR RNA004_HEK293_A RNA004_HEK293_B RNA004_UHRR_1 RNA004_UHRR_2 RNA004_blood_normal RNA004_blood_IVT \
    --outdir /home/awiercze/NanoComp_RNA002_RNA004_unaligned --raw --store --tsv_stats


# Run feature counts on gene level of all samples (using GENCODE gtf Release 43) 
bam_path_RNA002_RNA004=$( ls RNA00*_*GRCh38.bam )
featureCounts -a gencode.v43.annotation.gtf.gz -F 'GTF' -L -s 0 -T 64 -o RNA002_RNA004.GRCh38.featurecounts.tsv $RNA002_RNA004_aligned_bams

