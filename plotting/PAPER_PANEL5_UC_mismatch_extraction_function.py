import pysam 
import pandas as pd
import os
import sys
print ('argument list', sys.argv)
bam = sys.argv[1]

import pandas as pd
cols_mod_file = ["Chrom", #-> Chromosome.
"Start", #-> The starting coordinate.
"End", #-> The end coordinate.
"modID", #-> Unique ID identification.
"Score", #-> Score.
"Strand",# -> The direction of the strand containing this RNA modification site.
"modType", #-> The RNA modification type of this site.
"supportNum", #-> The number of datasets from which identified this RNA modification site.
"supportList", #-> The list of merged-datasets from which identified this RNA modification site.
"supportListSub", #-> The list of datasets from which identified this RNA modification site.
"pubList", #-> The list of pubmeds coorespond to above-mentioned datasets.
"cellList", #-> The list of cells or tissues in which this RNA modification is identified.
"seqTypeList", #-> The list of type of sequencing data of datasets.
"geneID", #-> The list of gene that the RNA modification site reside in.
"transcriptID", #-> The list of transcript that the RNA modification site reside in.
"geneName", #-> Gene name corresponding to the gene ID.
"geneType" ,#-> Gene biotype on which the RNA modification site distribute.
"Region", #-> Gene feature on which the RNA modification site distribute.
"Seq", #-> The sequence is 41-nt long that was extended by an additional 20 nt in both the 5′- and 3′-directions for the modification site.
"motifScore",# -> The 'Motif score' is alignment score to evaluate the accuracy of identified motif regions of m6A. The higher of motif score means a more accurate motif and a more reliable modification site. The range is from 0 to 5.
"conservedSitesList",# -> The list of sites that are conserved with this RNA modification site in other species.
"snoRNA_detailInfo", #-> Details of snoRNA that guide RNA modifications, eg. "U2:25|G|SCARNA2|C/D|snOPY" represents "snoRNA_guideSite", the base of this site, snoRNA name, snoRNA type and source of snoRNA.
"snoRNA_guideSite",# -> The position of the RNA modification site guided by snoRNAs on RNA, eg. "U2:25" means that this RNA modification site is in he 25th position of U2.
"snoRNA_nameList", #-> The list of snoRNAs that guide this RNA modification site.
"snoRNA_dataBase", #-> The source of information of this modification site guided by snoRNAs.
"writerLoc",# -> The genome coordinate of writer binding site. There may be multiple semicolon-separated values that correpond to the informations separated by semicolon in the 27-29 columns.
"writerID",# ->  The list of writer ID that is unique ID.
"writerNameList", #-> The name of writers that catalyze this RNA modification site.
"source"# -> The source ID (sample) of writer.
]

l = pd.read_csv("~/NEW_RNA004_PLOTS/human.hg38.Pseudo.result.col29.bed", sep = "\t", header = None, names = cols_mod_file)
l_mod = l[l["Chrom"].str.contains("chr")]



def extract_tavakoli_uc_mismatch(bam_file, tava_sites_in):
    p_all = pd.DataFrame()
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for i, row in tava_sites_in.iterrows():
        positions = [row[0], row[1], row[2]]
        n = 0
        C_mismatch = 0
        AUG_mismatch = 0
        for pileupcolumn in samfile.pileup(positions[0], start = positions[1], end = positions[2], truncate = True, max_depth = 10000000, min_base_quality=13): #defaul min_base_qualiy = 13
            out_df_mod = {"mismatch_C" : 0, "mismatch_AUG" : 0, "canonical" : 0, "deletion" : 0, "modified" : 0, "fail"  : 0}
            for pileupread in pileupcolumn.pileups:
                mod_pos = [m for m in pileupread.alignment.get_aligned_pairs(with_seq=True) if m[1] == pileupcolumn.pos]
                i = pileupread.alignment
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                if i.is_unmapped or i.is_supplementary or i.is_secondary:
                    continue
                if len(mod_pos) > 0:
                    if mod_pos[0][0] == None:
                        continue
                    n = n+1
                    if i.is_reverse:
                        if mod_pos[0][2] == "a":
                            called_base = i.query_sequence[mod_pos[0][0]]
                            if called_base == "G":
                                C_mismatch += 1
                            else:
                                AUG_mismatch += 1
                    else:
                        if mod_pos[0][2] == "t":
                            called_base = i.query_sequence[mod_pos[0][0]]
                            if called_base == "C":
                                C_mismatch += 1
                            else:
                                AUG_mismatch += 1
            new_df = pd.DataFrame({"sample" : os.path.basename(bam_file), "chrom" : positions[0], "chromStart" : positions[1], "n_reads" : n, "C_mismatch" : C_mismatch, "AUG_mismatch" : AUG_mismatch}, index = [0])
            p_all = pd.concat([p_all, new_df], axis = 0)
    return(p_all)



print(bam)
sample_name = os.path.basename(bam).replace(".bam", "")
p = extract_tavakoli_uc_mismatch(bam, tava_sites_in=l_mod)
p.to_csv("/home/awiercze/NEW_RNA004_PLOTS/RMBase_pseu_sites_UC_mismatch_"+ sample_name + ".csv", index = False)
