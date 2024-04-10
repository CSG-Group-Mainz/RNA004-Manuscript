#_____________________________________________________________________________________________
# RNA002 vs RNA004 SCRIPT: SNP CALLING OF VECTOR SEQUENCES OVER TIME 

# Load tools
ml pod5
ml dorado
ml samtools
ml minimap2
ml bcftools

# Functions
function run_dorado_RNA004 {
	time dorado basecaller \
        ./rna004_130bps_sup@v3.0.1 "$1"_basecall/"$2"_"$1"_basecall.pod5 \
        --modified-bases m6A_DRACH \
        --estimate-poly-a \
        -r \
       --device "cuda:0,cuda:1,cuda:2,cuda:3" \
       | samtools view -bh - \
       | samtools sort  -@ 32 - -o "$1"_basecall/"$2"_"$1"_basecall.unaligned.bam
}
export -f run_dorado_RNA004

function run_dorad_RNA002 {
	time dorado basecaller \
        ./rna002_70bps_hac@v3 "$1"_basecall/"$2"_"$1"_basecall.pod5 \
        --estimate-poly-a \
        -r \
        --device "cuda:0,cuda:1,cuda:2,cuda:3" \
       | samtools view -bh - \
       | samtools sort  -@ 32 - -o "$1"_basecall/"$2"_"$1"_basecall.unaligned.bam
}
export -f run_dorad_RNA002

function mapping {
	samtools fastq -T '*' "$1"_basecall/"$2"_"$1"_basecall.unaligned.bam \
	| minimap2 -y --MD -ax splice -uf -k14 -t 96 $3 - \
	| samtools view -bh \
	| samtools sort -@32 - -o "$1"_basecall/"$2"_"$1"_basecall.$4.bam
	samtools index "$1"_basecall/"$2"_"$1"_basecall.$4.bam
}
export -f mapping


function snp_call {
    n=$(expr $1 \* 4000)
	bcftools mpileup -r EGFP:115,mCherry:565 -a AD --max-depth $n -X ont -Ou -f $3 /raid/awiercze_analysis/RNA002_RNA004_paper/"$1"_basecall/"$2"_"$1"_basecall.$4.bam \
    | bcftools call -P 0.01 -m -A -Ov > /raid/awiercze_analysis/RNA002_RNA004_paper/SNP_CALL_OUTPUT/"$2"_"$1"_basecall.P1.$4.vcf
}
export -f snp_call

function create_input {
    echo "DONE"
    n=$(expr $1 \* 4000)
    echo $n
    mkdir -p "$1"_basecall
    tail -n+2 "$2"_seq_sum.txt | head -n $n | cut -f1 > "$1"_basecall/"$2"_"$1"_read_ids.txt 
}
export -f create_input

function filter_pod5 {
    pod5 filter $2 -t 96 -f --ids "$1"_basecall/"$3"_"$1"_read_ids.txt --output "$1"_basecall/"$3"_"$1"_basecall.pod5
}
export -f filter_pod5


# 0) Paths of pod5 files of HEK293+Vector Samples
RNA004_A="pod5_convert_A"
seq_sum_RNA004_A="DRS_004_A/20230727_1229_3A_PAQ16625_f4ea0e70/sequencing_summary_PAQ16625_f4ea0e70_dad6a4d3.txt"

RNA004_B="pod5_convert_B"
seq_sum_RNA004_B="DRS_004_B/20230727_1250_3B_PAQ68817_dccb2c17/sequencing_summary_PAQ68817_dccb2c17_548ac9d1.txt"

RNA002_A="DRS_002_A/20230714_1459_2A_PAG68216_8f015b5f/pod5/"
seq_sum_RNA002_A="DRS_002_A/20230714_1459_2A_PAG68216_8f015b5f/sequencing_summary_PAG68216_8f015b5f_2c064fb0.txt"

RNA002_B="DRS_002_B/20230714_1506_2B_PAG68254_d732e056/pod5/"
seq_sum_RNA002_B="DRS_002_B/20230714_1506_2B_PAG68254_d732e056/sequencing_summary_PAG68254_d732e056_88af490f.txt"


# 1) Sort sequencing summary file 
python3.8 sort_sequencing_summary.py $seq_sum_RNA004_A RNA004_A_seq_sum.txt
python3.8 sort_sequencing_summary.py $seq_sum_RNA004_B RNA004_B_seq_sum.txt
python3.8 sort_sequencing_summary.py $seq_sum_RNA002_A RNA002_A_seq_sum.txt
python3.8 sort_sequencing_summary.py $seq_sum_RNA002_B RNA002_B_seq_sum.txt

# 2) Create file directories with 10, 20, ... 100 files (=4000*x reads)
parallel -v -u --env create_input -j 10 create_input ::: {10..100..10} ::: "RNA004_A"
parallel -v -u --env create_input -j 10 create_input ::: {10..100..10} ::: "RNA004_B"
parallel -v -u --env create_input -j 10 create_input ::: {10..100..10} ::: "RNA002_A"
parallel -v -u --env create_input -j 10 create_input ::: {10..100..10} ::: "RNA002_B"

parallel -v -u --env create_input -j 1 create_input ::: 1 ::: "RNA004_A"
parallel -v -u --env create_input -j 1 create_input ::: 1 ::: "RNA004_B"
parallel -v -u --env create_input -j 1 create_input ::: 1 ::: "RNA002_A"
parallel -v -u --env create_input -j 1 create_input ::: 1 ::: "RNA002_B"


# 3) Create pod5 file of first 100 "files" to subset smaller pod5 files
i=100
pod5 filter $RNA004_A -t 96 -f --ids "$i"_basecall/RNA004_A_"$i"_read_ids.txt --output "$i"_basecall/RNA004_A_"$i"_basecall.pod5
pod5 filter $RNA004_B -t 96 -f --ids "$i"_basecall/RNA004_B_"$i"_read_ids.txt --output "$i"_basecall/RNA004_B_"$i"_basecall.pod5
pod5 filter $RNA002_A -t 96 -f --ids "$i"_basecall/RNA002_A_"$i"_read_ids.txt --output "$i"_basecall/RNA002_A_"$i"_basecall.pod5
pod5 filter $RNA002_B -t 96 -f --ids "$i"_basecall/RNA002_B_"$i"_read_ids.txt --output "$i"_basecall/RNA002_B_"$i"_basecall.pod5

# 3.1) Paths to pod5 file containing 100*4000 reads
pod5_100_RNA004_A="100_basecall/RNA004_A_100_basecall.pod5"
pod5_100_RNA004_B="100_basecall/RNA004_B_100_basecall.pod5"
pod5_100_RNA002_A="100_basecall/RNA002_A_100_basecall.pod5"
pod5_100_RNA002_B="100_basecall/RNA002_B_100_basecall.pod5"

# 4) Create pod5 files for 10 to 90 "files" by subsetting the 100 pod5 file
for i in {10..90..10}; do
    echo $i
    filter_pod5 $i $pod5_100_RNA004_A "RNA004_A"
    filter_pod5 $i $pod5_100_RNA004_B "RNA004_B"
    filter_pod5 $i $pod5_100_RNA002_A "RNA002_A"
    filter_pod5 $i $pod5_100_RNA002_B "RNA002_B"
done

# 4.1) 
filter_pod5 1 $pod5_100_RNA004_A "RNA004_A"
filter_pod5 1 $pod5_100_RNA004_B "RNA004_B"
filter_pod5 1 $pod5_100_RNA002_A "RNA002_A"
filter_pod5 1 $pod5_100_RNA002_B "RNA002_B"

# 5) Start basecalling, mapping and snp calling for PSU detection in vector sequences
for i in {10..100..10}; do
    echo $i
    # Basecall filtered reads
    run_dorado_RNA004 $i "RNA004_A" 
    run_dorado_RNA004 $i "RNA004_B"
    run_dorad_RNA002 $i "RNA002_A"
    run_dorad_RNA002 $i "RNA002_B"    

    # Map filtered reads
    mapping $i "RNA004_A" reference_vector.fa vector
    mapping $i "RNA004_B" reference_vector.fa vector
    mapping $i "RNA002_A" reference_vector.fa vector
    mapping $i "RNA002_B" reference_vector.fa vector


    # Call SNPs for modified region from mapped reads
    snp_call $i "RNA004_A" reference_vector.fa vector
    snp_call $i "RNA004_B" reference_vector.fa vector
    snp_call $i "RNA002_A" reference_vector.fa vector
    snp_call $i "RNA002_B" reference_vector.fa vector

done

# 5.1) Additionally, add results for 4000 reads
run_dorado_RNA004 1 "RNA004_A"
run_dorado_RNA004 1 "RNA004_B"
run_dorad_RNA002 1 "RNA002_A"
run_dorad_RNA002 1 "RNA002_B"    

mapping 1 "RNA004_A" reference_vector.fa vector
mapping 1 "RNA004_B" reference_vector.fa vector
mapping 1 "RNA002_A" reference_vector.fa vector
mapping 1 "RNA002_B" reference_vector.fa vector

snp_call 1 "RNA004_A" reference_vector.fa vector
snp_call 1 "RNA004_B" reference_vector.fa vector
snp_call 1 "RNA002_A" reference_vector.fa vector
snp_call 1 "RNA002_B" reference_vector.fa vector

## Plots of SNP calling results were built in R 

#_______________________________________________________________________________________________________________________
# GLOBAL BASECALLING ERROR PATTERN

# 1) Map dorado basecalled reads onto the EGFP and mCherry reference sequences using samtools and minimap2 and filter for primary reads only (-F3844)
for sample in *HEK293*unaligned.bam;do      
    echo  ${sample/unaligned.bam/vector.bam}
    samtools index -@64 $sample
	samtools fastq -T '*' $sample | minimap2 --MD -ax splice -uf -k14 -t 96 reporter_sequences.fa - | samtools view -hbS -F 3844 | samtools sort -@32 - -o ${sample/unaligned.bam/vector.bam}
	samtools index  ${sample/unaligned.bam/vector.bam}
done

# 2) Extract samtools mpileup
samtools mpileup -f reference_vector.fa -x --no-output-ins --no-output-del --no-output-ends -o RNA002_HEK293.mpileup RNA002_A_basecall.vector.bam  RNA002_B_basecall.vector.bam
samtools mpileup -f reference_vector.fa -x --no-output-ins --no-output-del --no-output-ends -o RNA004_HEK293.mpileup RNA004_A_basecall.vector.bam  RNA004_B_basecall.vector.bam
