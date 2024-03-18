# two arguments, a bam file and the tag to extract
BAM=$1
TAG=pt
OUT_PATH=$2

# write a tsv with columns read_id and tag value
echo -e "read_id\t$TAG"
samtools view "$BAM" | grep "$TAG:." | perl -pe 's/(^.+?)\t.*'$TAG':.:(.+?)\t.*/$1\t$2/g' > $OUT_PATH/poly_A_tail_est.tsv

#only keep length of unique reads
sort -s -k1,1 $OUT_PATH/poly_A_tail_est.tsv | uniq > $OUT_PATH/poly_A_tail_est_uniq.tsv
