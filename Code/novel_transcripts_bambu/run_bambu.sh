#!/bin/bash
source /raid/vincentd_analysis/miniforge3/etc/profile.d/conda.sh
conda activate bambu
ml samtools 

DDIR=/mnt/ssd_share_01/harmonized_rna004_test
BASEDIR=/raid/vincentd_analysis
REFDIR=$BASEDIR/references/NCBI

mkdir -p $BASEDIR/data
# echo "Merging 1 and 2 RNA004"
# samtools merge -@ 16 $BASEDIR/data/1_2_RNA004_sup.GRCh38.NCBIRefseq.bam $DDIR/1_RNA004_all_sup.GRCh38.NCBIRefseq.bam $DDIR/2_RNA004_all_sup.GRCh38.NCBIRefseq.bam
# echo "Merging A and B RNA002"
# samtools merge -f -@ 16 $BASEDIR/data/A_B_RNA002_hac.GRCh38.NCBIRefseq.bam $DDIR/A_RNA002_all_hac.GRCh38.NCBIRefseq.bam $DDIR/B_RNA002_all_hac.GRCh38.NCBIRefseq.bam
# echo "Merging A and B RNA004"
# samtools merge -f -@ 16 $BASEDIR/data/A_B_RNA004_sup.GRCh38.NCBIRefseq.bam $DDIR/A_RNA004_all_sup.GRCh38.NCBIRefseq.bam $DDIR/B_RNA004_all_sup.GRCh38.NCBIRefseq.bam

function run_bambu {
    BAMFILES=$1
    OUTDIR=$2

    mkdir -p $OUTDIR
    Rscript $BASEDIR/script_bambu.R \
        $BAMFILES \
        $REFDIR/GCF_000001405.40_GRCh38.p14_genomic.fna \
        $REFDIR/GCF_000001405.40_GRCh38.p14_genomic.gtf \
        $OUTDIR
}

run_bambu $DDIR/1_RNA002_all_hac.GRCh38.NCBIRefseq.bam,$BASEDIR/data/1_2_RNA004_sup.GRCh38.NCBIRefseq.bam $BASEDIR/bambu_results/1_2_RNA002_RNA004
run_bambu $BASEDIR/data/A_B_RNA002_hac.GRCh38.NCBIRefseq.bam,$BASEDIR/data/A_B_RNA004_sup.GRCh38.NCBIRefseq.bam $BASEDIR/bambu_results/A_B_RNA002_RNA004

conda deactivate
