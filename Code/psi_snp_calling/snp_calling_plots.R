################################################################################
####                                                                        ####
####    ANALYSIS OF SNP CALL OUTPUT FOR DIFFERENT SEQUENCING TIME POINTs    ####
####                                                                        ####
################################################################################


#_______________________________________________________________________________
# LOAD LIBRARIES
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

#_______________________________________________________________________________
# VCF INPUT FILES
# List all vcfs of different batches for RNA002 and RNA004 runs of Sample A and B
vcf_files = list.files("SNP_CALL_OUTPUT", ".vcf", full.names = T)

# Add names for each file
names(vcf_files) = gsub("_basecall.*", "", basename(vcf_files))

# Read all vcf files and create one data.frame of all samples using rbindlist
vcf_list = lapply(vcf_files, read_delim, delim = "\t", escape_double = FALSE, 
                                 col_names = c("CHROM",	"POS"	,"ID",	"REF",	"ALT","QUAL",	"FILTER",	"INFO",	"FORMAT",	"SAMPLE"), 
           comment = "#", trim_ws = TRUE)
vcf_df = data.table::rbindlist(vcf_list, idcol = "sample_id")  

# Extract AD column that contains the read depth for reference and alternate bases
vcf_df = separate(vcf_df, SAMPLE, into =  c("GT", "PL", "AD"), sep = ":")
vcf_df = separate(vcf_df, sample_id, into = c("FC", "Sample", "Batch"), remove = F)
vcf_df = separate(vcf_df, AD, into = c("U_count", "C_count"),sep = ",", remove = F)

# Calculate % of alternate base based on read depth for EGFP and mCherry at SNP positions
vcf_df$base_freq_alt = as.numeric(vcf_df$C_count)/(as.numeric(vcf_df$U_count)+as.numeric(vcf_df$C_count))
vcf_df$Significant = ifelse(vcf_df$GT %in% c("0/1", "1/1"), TRUE, FALSE)

#_______________________________________________________________________________
# EXTRACT TIME INFO FROM SEQ SUM FILES
# see python script for more information
seq_sum_100 = list.files(".", "[A|B]_seq_sum.txt", full.names = T)
names(seq_sum_100) = gsub("_seq_sum.*", "", basename(seq_sum_100))
seq_sum_list = lapply(seq_sum_100, function(i) {
  y = read.table(i, sep = "\t", header = T)
  print(nrow(y))
  out_df = data.frame()
  for (j in c(1, seq(10, 100, by = 10))){
    out_tmp = data_frame(batch = j, time = y[(j*4000), "end_time"])
    out_df = rbind(out_df, out_tmp)
  }
  out_df
})

seq_sum_df = data.table::rbindlist(seq_sum_list, idcol = "Sample_comp")
seq_sum_df$batch = as.character(seq_sum_df$batch)
seq_sum_df$sample_id = paste0(seq_sum_df$Sample_comp, "_", seq_sum_df$batch)

vcf_df = vcf_df %>% left_join(seq_sum_df, by = c("sample_id" = "sample_id", "Batch" = "batch"))

vcf_df$sample_id = gsub("_A", "_HEK293_A", vcf_df$sample_id )
vcf_df$sample_id = gsub("_B", "_HEK293_B", vcf_df$sample_id )
vcf_df$condition = gsub("_[A|B].*$", "", vcf_df$sample_id)
plot_df = unique(vcf_df[, c("sample_id", "Sample", "CHROM", "Batch", "base_freq_alt", "condition", "time", "Significant")])

#_______________________________________________________________________________
# VISUALIZE RESULTS FOR EACH BATCH AND SAMPLE
plot_df$sample_id_flat = gsub("_[1-9].*", "", plot_df$sample_id)
plot_df$colors = ifelse(plot_df$Significant, ifelse(plot_df$CHROM == "mCherry", "mCherry", "EGFP"), "not_signficant")

g = ggplot(plot_df, aes(x = time, 
                        y = base_freq_alt*100, 
                        color = CHROM, linetype = Significant, shape = Significant)) +
  
  geom_line() +
  geom_point(fill = "white") +
  facet_wrap(.~Sample + condition, ncol = length(unique(vcf_df$FC)), scales = "free_x") +
  theme_bw()+
  ylab("% bases U>C mismatch") +
  xlab("time [s]") +
  ylim(c(0,100)) +
  scale_color_manual("Vector", values = c("EGFP" = "#d90429", "mCherry" = "#8d99ae"))+ 
  scale_linetype_manual("SNP call", values= c("dashed", "solid")) +
  scale_shape_manual("SNP call", values = c(21,16)) +
  theme(legend.position = "bottom", axis.ticks.x = element_blank(),panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill="white"))

ggsave("SNPCalling_output_FINAL.svg", g, width = 9, height = 6, dpi = 300)
ggsave("SNPCalling_output_FINAL.png", g, width = 9, height = 6, dpi = 300)
ggsave("SNPCalling_output_FINAL.pdf", g, width = 9, height = 6, dpi = 300)




