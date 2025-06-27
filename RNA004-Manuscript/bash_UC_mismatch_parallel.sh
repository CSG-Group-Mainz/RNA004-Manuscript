#!/bin/bash 

# RUN U>C mismatch extraction in parallel 

bam_files=$( ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/DORADO_0_7_2_folder/*.0.7.2.GRCh38.bam)
bam_files2=$( ls /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S*/RNA004_S*.0.7.2.GRCh38.bam)
bam_files3=$( ls /raid/awiercze_analysis/RNA00*_UHRR_basecall.GRCh38.bam /raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S*/RNA004_S*.0.7.2.GRCh38.bam)



function extract_UC_mismatch {
    echo $1
    python3.8 run_UC_mismatch_extraction.py $1
}
export -f extract_UC_mismatch

#parallel -u --env extract_UC_mismatch -j 27 extract_UC_mismatch ::: $bam_files
#parallel -u --env extract_UC_mismatch -j 27 extract_UC_mismatch ::: $bam_files2
parallel -u --env extract_UC_mismatch -j 27 extract_UC_mismatch ::: $bam_files3
