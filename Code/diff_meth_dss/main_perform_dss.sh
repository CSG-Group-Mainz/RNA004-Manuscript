#!/bin/bash
source /home/vincent/miniforge3/etc/profile.d/conda.sh
BASEDIR=/home/vincent/projects/charlotte_rna2_4/differential_methylation_analysis

# modkit pileup files must be placed in $BASEDIR/modkit_pileup
echo "Prepping modkit pileup files"
for i in $BASEDIR/modkit_pileup/*.GRCh38.bam.bed;
do
   python $BASEDIR/script_prep_pileup_for_dss.py --path $i
done

echo "Running DSS - forward strand positions only"
mkdir -p $BASEDIR/dss_output
Rscript $BASEDIR/script_run_dss_tworep.R \
    $BASEDIR/modkit_pileup/RNA004_A_basecall.GRCh38.bam_forw_prepped.bed \
    $BASEDIR/modkit_pileup/RNA004_B_basecall.GRCh38.bam_forw_prepped.bed \
    $BASEDIR/modkit_pileup/RNA004_UHRR_1_basecall.GRCh38.bam_forw_prepped.bed \
    $BASEDIR/modkit_pileup/RNA004_UHRR_2_basecall.GRCh38.bam_forw_prepped.bed \
    "HEK293T" "UHRR" "forward" \
    $BASEDIR/dss_output

echo "Running DSS - reverse strand positions only"
Rscript $BASEDIR/script_run_dss_tworep.R \
    $BASEDIR/modkit_pileup/RNA004_A_basecall.GRCh38.bam_rev_prepped.bed \
    $BASEDIR/modkit_pileup/RNA004_B_basecall.GRCh38.bam_rev_prepped.bed \
    $BASEDIR/modkit_pileup/RNA004_UHRR_1_basecall.GRCh38.bam_rev_prepped.bed \
    $BASEDIR/modkit_pileup/RNA004_UHRR_2_basecall.GRCh38.bam_rev_prepped.bed \
    "HEK293T" "UHRR" "reverse" \
    $BASEDIR/dss_output

echo "Transforming DSS output files to bed-like format"
python $BASEDIR/script_dss_output_to_bed.py --path $BASEDIR/dss_output/forward_HEK293T_UHRR_dml.tsv
python $BASEDIR/script_dss_output_to_bed.py --path $BASEDIR/dss_output/reverse_HEK293T_UHRR_dml.tsv

echo "Transforming GFF to bed-like format"
# downloaded into $BASEDIR/gff and unzipped beforehand
python $BASEDIR/script_gff_to_bed.py --path /home/vincent/projects/charlotte_rna2_4/differential_methylation_analysis/gff/gencode.v45.primary_assembly.basic.annotation.gff3

echo "Intersecting DSS output with GFF"
mkdir -p $BASEDIR/dss_output_annotated
bedtools intersect -wao -a $BASEDIR/dss_output/forward_HEK293T_UHRR_dml.bed -b $BASEDIR/gff/gencode.v45.primary_assembly.basic.annotation_forward.tsv > $BASEDIR/dss_output_annotated/forward_HEK293T_UHRR_dml.tsv
bedtools intersect -wao -a $BASEDIR/dss_output/reverse_HEK293T_UHRR_dml.bed -b $BASEDIR/gff/gencode.v45.primary_assembly.basic.annotation_reverse.tsv > $BASEDIR/dss_output_annotated/reverse_HEK293T_UHRR_dml.tsv

echo "Processing outputs and creating plots"

mkdir -p $BASEDIR/plots

conda activate rna4_2
python $BASEDIR/script_create_plots.py \
    --path_f $BASEDIR/dss_output_annotated/forward_HEK293T_UHRR_dml.tsv \
    --path_r $BASEDIR/dss_output_annotated/reverse_HEK293T_UHRR_dml.tsv \
    --path_f_annot $BASEDIR/dss_output_annotated/forward_HEK293T_UHRR_dml.tsv \
    --path_r_annot $BASEDIR/dss_output_annotated/reverse_HEK293T_UHRR_dml.tsv \
    --outdir $BASEDIR/plots

conda deactivate