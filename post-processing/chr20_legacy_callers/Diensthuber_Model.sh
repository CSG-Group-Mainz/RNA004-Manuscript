#this scripts lists the commands needed to rebasecall the chr20 data with the models from the diensthuber publication, that were downloaded as described there.
#it is furthermore assumed at this point that the data has been subset to chr20 and that pod5 has been converted to fast5, as described in subset_chr20.sh. because running all the models on all the data is ressource and time-prohibitive.

#load guppy
ml guppy/6.2.1
ml samtools/1.16.1
ml minimap2/2.26

echo "run guppy w/ diensthuber models"
#call the ivt model for ivt data
guppy_basecaller -i ./chr20_rna002_blood_IVT_fast5/ -s ./Diensthuber_SUP_IVT_guppy_out_blood -c /home/chhewel/rna_r9.4.1_70bps_ivt_sup.cfg -m /home/chhewel/template_rna_r9.4.1_70bps_ivt_sup.jsn -x cuda:all --bam_out


#call the sup model for our hek293t dat,  blood data and uhrr
guppy_basecaller -i ./chr20_rna002_lemke_b_fast5/ -s ./Diensthuber_SUP_guppy_out_B -c /home/chhewel/rna_r9.4.1_70bps_sup.cfg -m /home/chhewel/template_rna_r9.4.1_70bps_sup.jsn -x cuda:all --bam_out
guppy_basecaller -i ./chr20_rna002_lemke_a_fast5/ -s ./Diensthuber_SUP_guppy_out_A -c /home/chhewel/rna_r9.4.1_70bps_sup.cfg -m /home/chhewel/template_rna_r9.4.1_70bps_sup.jsn -x cuda:all --bam_out
guppy_basecaller -i ./chr20_rna002_blood_fast5/ -s ./Diensthuber_SUP_guppy_out_blood -c /home/chhewel/rna_r9.4.1_70bps_sup.cfg -m /home/chhewel/template_rna_r9.4.1_70bps_sup.jsn -x cuda:all --bam_out
guppy_basecaller -i ./chr20_rna002_uhrr_fast5/ -s ./Diensthuber_SUP_guppy_out_UHRR -c /home/chhewel/rna_r9.4.1_70bps_sup.cfg -m /home/chhewel/template_rna_r9.4.1_70bps_sup.jsn -x cuda:all --bam_out



#align against genome
echo "step1 map the stuff"
function map2GRCh38 {
	small=$( echo $1 | cut -f2 -d "/" )
	echo $small
	          minimap2  -y --MD -ax splice -uf -k14 -t 32 GRCh38.primary_assembly.genome.fa.gz $1/*fastq | samtools sort  -@ 16 - -o ${small}.GRCh38.bam
	
	samtools index -@8 ${small}.GRCh38.bam

	}
 
export -f map2GRCh38
bam=$( ls *bam)
RNA004_unaligned_bams=$( ls -d ./Diensthuber_*A/pass )
parallel -v -u --env map2GRCh38 -j 4 map2GRCh38 ::: $RNA004_unaligned_bams

#run pomoxis stats_from_bam and nanocomp to extract basic QC information
echo "now pomoxis"
function get_pomoxis {
	 stats_from_bam  -o ${1/.bam/.summary-txt} -s ${1/.bam/.pomoxis.txt} -t 16 $1
	         }

	 export -f get_pomoxis

	 parallel -v -u --env get_pomoxis -j 4 get_pomoxis ::: $bam
	 
#run Nanocomp
echo "now nanocomp"
bam_list=$(ls *.bam  )
NanoComp --bam $bam_list -t 40  --outdir ./Diensthuber_NanoComp_3 --raw --store --tsv_stats

