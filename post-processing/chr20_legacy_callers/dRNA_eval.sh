#the purpose of this script is to run the dRNA eval code from here https://github.com/KleistLab/nanopore_dRNAseq
#for this we take 1) a subsection of the data that aligns to chr20 and 2) re-align the file before processing, because they request specific options such as SAM format, discarding secondary hits, which would not be feasible for the main bams used elsewhere
ml samtools
ml minimap2
#get list of input files chr 20 bams 
input=$( ls *20.bam )
#okay transfer file format into input data for eval.py and subsequently use eval.py
function do_eval {
	echo $1
	#remap according to specs in kleistlab github
samtools fastq -@ 32  -T'*' $1 | minimap2 -t 32 -a -secondary=no --eqx --sam-hit-only GRCh38.primary_assembly.genome.fa - >  ${1/.bam/.sam}
#record which sample was processed
echo -e "sample\t$1" > eval.stats.$1.txt
#run eval and store output in stats file
 python ~/NanoPore_dRNASeq/eval.py GRCh38.primary_assembly.genome.fa  ${1/.bam/.sam} >> eval.stats.$1.txt
 #delete samfile again to minimize storage volume used
 rm ${1/.bam/.sam}

}
export -f do_eval
parallel -j4 -u do_eval ::: $input
#extract group specific stats and merge them into one file for easier handling
for i in $( ls eval.stats* ) ; do echo $i ; name=$( echo $i | sed 's/eval.stats.//; s/.0.7.2.GRCh38.chr20.bam.txt//; s/_basecall//' ) ;  grep acc $i | sed "s/The overall accuracy on this dataset is/accuracy/; s/ (/\t/; s/ i/\ti/;s/ d/\td/; s/)/\t$name/" >> merged_accuracy_eval.stats.txt; done
for i in $( ls eval.stats* ) ; do echo $i ; name=$( echo $i | sed 's/eval.stats.//; s/.0.7.2.GRCh38.chr20.bam.txt//; s/_basecall//' ) ;  grep "sequencing errors" $i | sed "s/The sequencing errors are //; s/ and/\t/;  s/homopolymers/homopolymers\t$name/" >> merged_sequencing_errors_eval.stats.txt; done

#make everything tsv for easier reading in sed -i does in place edit of files
for i in merged*txt; do echo $i; sed -i 's/ /\t/g' $i ; done
