ml pod5
ml dorado
ml samtools
ml minimap2

export FILEPATH_IVT=$1
export FILEPATH_DIR=$2
export VAR=$3
export REF_GENOME="/raid/fehof_analysis/ref_genome/GRCh38.primary_assembly.genome.fa"
export REF_ANNOT="/raid/fehof_analysis/ref_genome/gencode.v43.basic.annotation.gff3"
export MODKIT="/raid/fehof_analysis/software/modkit/modkit"
export FEATURE_COUNTS="/mnt/ssd_share_01/harmonized_rna004_test/featureCounts"
export GENE_TARGETS="/raid/fehof_analysis/parallel_test/Direct_RNA004/gene_targets.txt"

#basecall the IVT data with m6A basecalling model + poly(A) tail estimation
	time dorado basecaller \
        /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 $FILEPATH_IVT \
        --modified-bases m6A_DRACH \
        --estimate-poly-a \
        --device "cuda:0,cuda:1,cuda:2,cuda:3" \
        | samtools view -bh - \
        | samtools sort  -@ 32 - -o $FILEPATH_IVT/$VAR.basecall.unaligned.bam

#basecall the Direct data with m6A basecalling model + poly(A) tail estimation
	time dorado basecaller \
        /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 $FILEPATH_DIR \
        --modified-bases m6A_DRACH \
        --estimate-poly-a \
        --device "cuda:0,cuda:1,cuda:2,cuda:3" \
        | samtools view -bh - \
        | samtools sort  -@ 32 - -o $FILEPATH_DIR/$VAR.basecall.unaligned.bam

#alignment of basecalled data to hg38 reference
function mapping {
	#echo $1/$VAR.basecall.unaligned.bam
	samtools fastq -T '*' $1/$VAR.basecall.unaligned.bam \
	| minimap2 -y --MD -ax splice -uf -k14 -t 32 $REF_GENOME - \
	| samtools view -bh - \
	| samtools sort -@32 - -o $1/$VAR.hg38.aligned.bam
	samtools index $1/$VAR.hg38.aligned.bam
}

export -f mapping

parallel -v -u --env mapping -j 2 mapping ::: $FILEPATH_IVT $FILEPATH_DIR

#create bed file for loreme and visualize in a metagene plot for m6A peaks
function metagene_plot {
	$MODKIT pileup -t 32 $1/$VAR.hg38.aligned.bam $1/$VAR.GRCh38_aligned.bed --log-filepath $1/pileup-log
#	loreme mean --total $1/$VAR.GRCh38_aligned.bed
	loreme plot-genes --smooth --title "Methylation plot for $VAR files" $REF_ANNOT $1/$VAR.GRCh38_aligned.bed $1/$VAR.genebody_methylation.RNA004_GRCh38_aligned.png
}

export -f metagene_plot

parallel -v -u --env metagene_plot -j 2 metagene_plot ::: $FILEPATH_IVT $FILEPATH_DIR

function featureCounts {
	$FEATURE_COUNTS -T 32 --verbose -F 'GTF' -L -s 0 -a $REF_ANNOT -o $1/$VAR.featureCounts.RNA004.hg38.csv $1/$VAR.hg38.aligned.bam
	cat $1/$VAR.featureCounts.RNA004.hg38.csv | grep -f $GENE_TARGETS > $1/$VAR.featureCounts.RNA004.hg38_filtered.csv
}

export -f featureCounts

parallel -v -u --env featureCounts -j 2 featureCounts ::: $FILEPATH_IVT $FILEPATH_DIR

function poly_a_estimation {

        bash /raid/fehof_analysis/scripts/extract_poly_a_tag.sh $1/$VAR.hg38.aligned.bam $1
}

export -f poly_a_estimation

parallel -v -u --env poly_a_estimation -j 2 poly_a_estimation ::: $FILEPATH_IVT $FILEPATH_DIR
