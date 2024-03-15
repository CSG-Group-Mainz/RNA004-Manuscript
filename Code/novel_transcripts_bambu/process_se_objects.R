##### 1_2 #####
library(bambu)

# load object
se = readRDS("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/1_2_RNA002_RNA004/se_object.rds")

se.novel = se[mcols(se)$novelTranscript,]
writeBambuOutput(se.novel, "/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/1_2_RNA002_RNA004/1_2_novel")

se.fullen = se[(assays(se)$fullLengthCounts >= 1)[,1] | (assays(se)$fullLengthCounts >= 1)[,2],]
writeBambuOutput(se.fullen, "/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/1_2_RNA002_RNA004/1_2_fulllen")


##### A_B #####

se = readRDS("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/A_B_RNA002_RNA004/se_object.rds")

se.novel = se[mcols(se)$novelTranscript,]
writeBambuOutput(se.novel, "/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/A_B_RNA002_RNA004/A_B_novel")

se.fullen = se[(assays(se)$fullLengthCounts >= 1)[,1] | (assays(se)$fullLengthCounts >= 1)[,2],]
writeBambuOutput(se.fullen, "/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/A_B_RNA002_RNA004/A_B_fulllen")
