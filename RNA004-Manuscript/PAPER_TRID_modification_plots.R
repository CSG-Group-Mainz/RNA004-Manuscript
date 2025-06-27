#________________________________________________________________________________#
#               VECTOR MODIFICATION PLOTS                                        #
#________________________________________________________________________________#
x = c("modified" = "#B39BC8", "unmodified" = "#d48a79")

#________________________________________________________________________________
# LOAD LIBRARIES
library(ggplot2)
library(data.table)
library(readr)
library(dplyr)

#________________________________________________________________________________
# SET WORKING DIRECTORY THAT CONTAINS OUTPUT OF JUPYTER NOTEBOOK 
setwd("/raid/awiercze_analysis/RESEARCH/REVISION_FIGURES")

#________________________________________________________________________________
# LOAD CSV FILES CONTAINING THE NUMBER OF READS WITH PSEU CALLS ON U OR C, U>C 
# MISMATCHES AND UNMODIFIED CALLS FOR THE POSITIONS OF INTEREST IN THE SYNTHETIC
# OLIGOS, REPORTERS AND PSMB2 GENE FROM HEK293T CELLS

# List files that starts with "Paper_pseU_calls_" in current dir
modis = list.files(".", pattern = "PAPER_pseU_calls_")
modis = modis[grep(".csv", modis)]
# Add Source labels
names(modis) = c("Synthetic_oligo","PSMB2", "mRNA_vector")
# Read files 
modi_list = lapply(modis, function(i) {
    x = fread(i)
    x = x[, c("condition", "NreadsU_calls", "ref_pos","Ref", "Contig", "NreadsC_calls")]
    x
})
# Combine dataframes
modi_df = rbindlist(modi_list, idcol = "Source")

# Rename Position Labels per Sample
modi_df[["Sample"]] = modi_df$Ref
modi_df[["Sample"]][modi_df$Ref == "EGFP" & modi_df$ref_pos == 73] = "Oligo1_EGFP_motif"
modi_df[["Sample"]][modi_df$Ref == "EGFP" & modi_df$ref_pos == 113] = "Oligo2_EGFP_motif"
modi_df[["Sample"]][modi_df$Ref == "mCherry" & modi_df$ref_pos == 73] = "Oligo2_mCherry_motif"
modi_df[["Sample"]][modi_df$Ref == "mCherry" & modi_df$ref_pos == 113] = "Oligo1_mCherry_motif"
modi_df[["Sample"]][modi_df$Ref == "EGFP" & modi_df$Contig == "Sample A"] = "SampleA_EGFP"
modi_df[["Sample"]][modi_df$Ref == "EGFP" & modi_df$Contig == "Sample B"] = "SampleB_EGFP"
modi_df[["Sample"]][modi_df$Ref == "mCherry" & modi_df$Contig == "Sample A"] = "SampleA_mCherry"
modi_df[["Sample"]][modi_df$Ref == "mCherry" & modi_df$Contig == "Sample B"] = "SampleB_mCherry"
modi_df[["Sample"]] = gsub("RNA004_", "", modi_df[["Sample"]])
modi_df[["Sample"]] = gsub("A$", "HEK293T_A", modi_df[["Sample"]])
modi_df[["Sample"]] = gsub("B$", "HEK293T_B", modi_df[["Sample"]])

# Convert columns of modification status, sample and source to factor levels
#modi_df$modified = factor(modi_df$condition, levels = rev(c("modified", "unmodified")))
modi_df$Sample = factor(modi_df$Sample, levels = c("SampleA_EGFP", "SampleA_mCherry", "SampleB_EGFP", "SampleB_mCherry", "HEK293T_A", "HEK293T_B", "blood", "blood_IVT", "UHRR_1", "UHRR_2", "UHRR_3", "Oligo1_EGFP_motif", "Oligo1_mCherry_motif", "Oligo2_mCherry_motif", "Oligo2_EGFP_motif"))
modi_df$Source = factor(modi_df$Source, levels = c("Synthetic_oligo","mRNA_vector", "PSMB2")) 

# Subtract the number of pseU calls on C sites from the U>C mismatch counts of U-bases pseU calls 
# and add a new condition containing C-based pseU calls
modi_df_modified = modi_df[modi_df$condition == "modified",]

modi_df_modified_tmp = modi_df_modified
modi_df_modified_tmp$condition = "modified_C_calls"
modi_df_modified_tmp$NreadsU_calls = modi_df_modified_tmp$NreadsC_calls
modi_df_modi_calls = rbind(modi_df_modified, modi_df_modified_tmp)
modi_df_modi_calls = modi_df_modi_calls[, c("Source", "condition", "NreadsU_calls", "ref_pos", "Ref", "Contig", "Sample")]


modi_df_C_mismatch = modi_df[modi_df$condition == "mismatch_C",]
tmp = merge(modi_df_modified[, c("Source", "condition", "NreadsC_calls", "ref_pos", "Ref", "Contig", "Sample")], 
      modi_df_C_mismatch[, c("Source", "condition", "NreadsU_calls", "ref_pos", "Ref", "Contig", "Sample")],
      by = c("Source", "ref_pos", "Ref", "Contig", "Sample"))
tmp$NreadsU_calls = tmp$NreadsU_calls - tmp$NreadsC_calls
tmp = tmp[, c("Source", "condition.y", "NreadsU_calls", "ref_pos", "Ref", "Contig", "Sample")]
colnames(tmp) = c("Source", "condition", "NreadsU_calls", "ref_pos", "Ref", "Contig", "Sample")
modi_df_modi_calls = rbind(modi_df_modi_calls, tmp)

# Label canonical bases and bases that did not reach the probability threshold for modification as "unmodified"
# and summarize per sample
modi_df_tmp_unmodified = modi_df[modi_df$condition %in% c("canonical", "fail"),]
modi_df_tmp_unmodified_f = modi_df_tmp_unmodified %>% group_by(Source,ref_pos, Ref, Contig, Sample) %>% 
    summarise(NreadsU_calls_sum = sum(NreadsU_calls)) %>%
    mutate(condition = "unmodified")
modi_df_tmp_unmodified_f = modi_df_tmp_unmodified_f[, c("Source", "condition", "NreadsU_calls_sum", "ref_pos", "Ref", "Contig", "Sample")]
colnames(modi_df_tmp_unmodified_f) = c("Source", "condition", "NreadsU_calls", "ref_pos", "Ref", "Contig", "Sample")

# Combine daraframes of modified calls and unmodified calls
modi_df_modi_calls = rbind(modi_df_modi_calls, as.data.frame(modi_df_tmp_unmodified_f))

# Calculate the frequency of each condition over all counts
modi_df_new_perc = modi_df_modi_calls %>% group_by(Source, ref_pos, Ref, Contig) %>% 
    mutate("%" = NreadsU_calls/sum(NreadsU_calls)*100)

# Define order of conditions in plot by adding levels 
modi_df_new_perc$condition = factor(modi_df_new_perc$condition, levels = rev(c("modified", "modified_C_calls", "mismatch_C", "unmodified")))
modi_df_new_perc$Sample = factor(modi_df_new_perc$Sample, levels = c("SampleA_EGFP", "SampleA_mCherry", "SampleB_EGFP", "SampleB_mCherry", "HEK293T_A", "HEK293T_B", "Oligo1_EGFP_motif", "Oligo2_EGFP_motif", "Oligo2_mCherry_motif", "Oligo1_mCherry_motif"))

# Define color palette
0b3954-098c9a-a0a0a0--8b5cf6
final_palette <- rev(c("#ee6c4d", "#F5A38F", "#8B5CF6", "lightgray"))

modi_df_new_perc_sub = modi_df_new_perc %>% filter(Sample %in% c("SampleA_EGFP", "SampleA_mCherry", "SampleB_EGFP", "SampleB_mCherry", "HEK293T_A", "HEK293T_B", "Oligo1_EGFP_motif", "Oligo2_EGFP_motif", "Oligo2_mCherry_motif", "Oligo1_mCherry_motif"))
#Plot frequencies per condition and sample
g = ggplot(modi_df_new_perc_sub, aes(x = Sample, y = `%`, fill = condition) )+
    geom_bar(stat = "identity", color = "black") +
    facet_grid(.~Source, scales = "free", space = "free") +
    scale_x_discrete(expand = c(0.1, 0.4)) +
    theme_classic() +
    scale_fill_manual("", values = final_palette) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    ylab("% pseU modified")
g
ggsave("/home/awiercze/PAPER_FIGURE_5F.pdf", g, width = 7, height = 5, dpi = 300)
ggsave("/home/awiercze/PAPER_FIGURE_5F.pdf", g, width = 7, height = 5, dpi = 300)

# Save table of read counts per condition and sample
write.table(modi_df_new_perc_sub, "/raid/awiercze_analysis/RESEARCH/REVISION_FIGURES/PAPER_SUPP_TABLE_S4.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#____________________________________________________________________________
# Define U>C mismatch and dorado pseU calls on U and C-sites as modified and some counts
modi_df_modi_calls$condition2 = ifelse(modi_df_modi_calls$condition == "unmodified", "unmodified", "modified")
modi_df_modi_calls_sub2 = modi_df_modi_calls %>% group_by(Source,ref_pos, Ref, Contig, Sample, condition2) %>% 
    summarise(NreadsU_calls_sum = sum(NreadsU_calls))

# Calculate frequency for unmodified or modified counts
modi_df_new2_perc = modi_df_modi_calls_sub2 %>% group_by(Source, ref_pos, Ref, Contig) %>% 
    mutate("%" = NreadsU_calls_sum/sum(NreadsU_calls_sum)*100)

# Define levels for plot order
modi_df_new2_perc$condition2 = factor(modi_df_new2_perc$condition2, levels = rev(c("modified", "unmodified")))
modi_df_new2_perc$Sample = factor(modi_df_new2_perc$Sample, levels = c("SampleA_EGFP", "SampleA_mCherry", "SampleB_EGFP", "SampleB_mCherry", "HEK293T_A", "HEK293T_B", "Oligo1_EGFP_motif", "Oligo2_EGFP_motif", "Oligo2_mCherry_motif", "Oligo1_mCherry_motif"))
modi_df_new2_perc$Source = factor(modi_df_new2_perc$Source, levels = c("Synthetic_oligo","mRNA_vector", "PSMB2")) 

modi_df_new2_perc_sub = modi_df_new2_perc %>% filter(Sample %in% c("SampleA_EGFP", "SampleA_mCherry", "SampleB_EGFP", "SampleB_mCherry", "HEK293T_A", "HEK293T_B", "Oligo1_EGFP_motif", "Oligo2_EGFP_motif", "Oligo2_mCherry_motif", "Oligo1_mCherry_motif"))
# Plot proportion of modified and unmodified reads per sample
g = ggplot(modi_df_new2_perc_sub, aes(x = Sample, y = `%`, fill = condition2) )+
    geom_bar(stat = "identity", color = "black") +
    facet_grid(.~Source, scales = "free", space = "free") +
    scale_x_discrete(expand = c(0.1, 0.4)) +
    theme_classic() +
    scale_fill_manual("", values = c("lightgray","#ee6c4d")) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    ylab("% pseU modified")
g

ggsave("/home/awiercze/PAPER_FIGURE_5G.pdf", g, width = 7, height = 5, dpi = 300)
ggsave("/home/awiercze/PAPER_FIGURE_5G.png", g, width = 7, height = 5, dpi = 300)


