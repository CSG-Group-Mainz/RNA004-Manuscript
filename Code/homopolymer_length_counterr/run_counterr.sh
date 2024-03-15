#!/bin/bash
source /raid/vincentd_analysis/miniforge3/etc/profile.d/conda.sh

DDIR=/mnt/ssd_share_01/harmonized_rna004_test
ODIR=/raid/vincentd_analysis/counterr_output
REF=/raid/vincentd_analysis/references/NCBI/GCF_000001405.40_GRCh38.p14_genomic.fna

conda activate counterr

function run_counterr {
	BAM=$1
	OUTDIR=$ODIR/$2

	echo "Starting counterr: $BAM | Writing to: $OUTDIR"

	mkdir -p $OUTDIR
	counterr -bam $BAM -genome $REF -output_dir $OUTDIR
}


run_counterr $DDIR/1_RNA002_all_hac.GRCh38.NCBIRefseq.bam 1_RNA002_GRCh38_NCBIRefseq
run_counterr $DDIR/1_RNA004_all_sup.GRCh38.NCBIRefseq.bam 1_RNA004_GRCh38_NCBIRefseq
run_counterr $DDIR/2_RNA004_all_sup.GRCh38.NCBIRefseq.bam 2_RNA004_GRCh38_NCBIRefseq

run_counterr $DDIR/A_RNA002_all_hac.GRCh38.NCBIRefseq.bam A_RNA002_GRCh38_NCBIRefseq
run_counterr $DDIR/A_RNA004_all_sup.GRCh38.NCBIRefseq.bam A_RNA004_GRCh38_NCBIRefseq

run_counterr $DDIR/B_RNA002_all_hac.GRCh38.NCBIRefseq.bam B_RNA002_GRCh38_NCBIRefseq
run_counterr $DDIR/B_RNA004_all_sup.GRCh38.NCBIRefseq.bam B_RNA004_GRCh38_NCBIRefseq

run_counterr $DDIR/Control_RNA002_hac.GRCh38.NCBIRefseq.bam Control_RNA002_GRCh38_NCBIRefseq

conda deactivate
