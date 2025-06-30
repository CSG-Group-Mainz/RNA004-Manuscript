
ml modkit/0.3.1

#the input for this step are the basecalled GRCh38 aligned bams with MM/ML tags for m6a and pseu from the script in main and the output are filtered mod bed files
#the purpose of this file is to extract modification stochimetry info in accordance with 1)  mod threshold and 2) coverage threshold; 
#the mod threshold was derived from the modifications histograms of IVT and normal samples; the coverage and % methylation filter was derived from the maximum of F1-score setting the GLORI data as ground truth; the DRACH filter gets applied to m6A because the all context model was used during the basecall step and not filtering to DRACH regions would negatively impact precision and recall for true sites on  m6A samples  

#step 0 gather samples 
bams_to_use=$( ls ./RNA004_Manuscript_reanalysis/*38.bam )

echo "start pileup"


#step 1 pileup modkit
function do_pileup {
               echo $1

                       base=$(basename $1)
		       echo $base
                       modkit pileup --threads 18 -r GRCh38.primary_assembly.genome.fasta --ignore a  --log-filepath ${base/.bam/.modkit.log}  --filter-threshold T:0.8 --mod-threshold 17802:0.98 $1 ${base/.bam/_pseu.r1.mod.bed}
			modkit pileup --threads 16 --motif DRACH 2 -r GRCh38.primary_assembly.genome.fasta --ignore 17802  --log-filepath ${base/.bam/.modkit.log}  --filter-threshold A:0.8 --mod-threshold a:0.98 $1 ${base/.bam/_m6A.r1.mod.bed}
 			
			modkit sample-probs --threads 8 $1 --hist -o hist_dir${base/.bam/}  --percentiles 0.1,0.8,0.85,0.9
#
                       }
               export -f do_pileup

               parallel -v do_pileup ::: $bams_to_use

#step2 filter to coverage and percent methylation threshold

mod_beds_to_filter=$( ls *.r1.mod.bed )
function do_filter {
	           echo $1
		   awk -F "\t" '{ if(($11 >= 10 ) && ($5 >= 10 )) { print } }' <$1 >${1/.mod.bed/.filter10_10.mod.bed}
		export -f do_filter
	parallel -v do_filter ::: $mod_beds_to_filter
#
