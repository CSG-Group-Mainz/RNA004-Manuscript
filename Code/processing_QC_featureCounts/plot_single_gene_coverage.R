################################################################################
####                                                                        ####
####                     SINGLE GENE BODY COVERAGE PLOTS                    ####
####                                                                        ####
################################################################################

#_______________________________________________________________________________
# LOAD LIBRARIES
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

#_______________________________________________________________________________
# DEFINE COLOR SCHEME
colors_input = c("RNA002_UHRR" = "#005f73", 
"RNA002_HEK293" =  "#0a9396",
"RNA002_blood" ="#94d2bd", 
"RNA004_blood" = "#ee9b00", 
"RNA004_HEK293" = "#ca6702", 
"RNA004_UHRR" = "#bb3e03")
theme_set(theme_classic())

#_______________________________________________________________________________
# DEFINE GENES OF INTEREST AND OUTPUT FOLDER FROM RSEQC
gene_list =c("XIST", "MT-ND3")

input_dir = "/raid/awiercze_analysis/RNA002_RNA004_paper"

# RUN PLOTTING FUNCTION
for (geneOfInterest in gene_list){
  #_______________________________________________________________________________
  # REPLACE - with _
  geneOfInterest = gsub("-", "_", geneOfInterest)

  #_______________________________________________________________________________
  # SOURCE R SCRIPTS FROM genebody_coverage.py FUNCTION OUTPUT FOR GENE OF INTEREST
  # (counts per percentiles are stored in .r files)
  samples_file = list.files(input_dir, pattern = paste0(geneOfInterest, ".geneBodyCoverage.r"), full.names = T)
  lapply(samples_file, source)
  print(samples_file)
  #_______________________________________________________________________________
  # Extract all loaded variables 
  file_list = ls()[grep(geneOfInterest, ls())]
  print(file_list)
  names(file_list) = gsub(paste0("_basecall.GRCh38.", geneOfInterest), "", file_list)
  names(file_list) = gsub("_A", "_HEK293_A", names(file_list))
  names(file_list) = gsub("_B", "_HEK293_B", names(file_list))

  geneBody_list = lapply(file_list, function(i) {
      p = get(i)
      q = data.frame(Percentiles = 1:100, PercVal = p)
      q
  })

  geneBodyCov = rbindlist(geneBody_list, idcol = "ID")


  #_______________________________________________________________________________
  # ADD CONDITIONS TO TABLE FOR COLOR AND SHAPE PLOTTING
  geneBodyCov$Condition = gsub("_[A|B|D|I|1|2].*", "", geneBodyCov$ID)
  geneBodyCov$Shape = 1
  geneBodyCov$Shape = ifelse(geneBodyCov$ID %in% c("RNA002_HEK293_B", "RNA004_HEK293_B", "RNA002_blood_IVT", "RNA004_blood_IVT", "RNA004_UHRR_2"), "2", "1")

  n_samples = length(unique(geneBodyCov$ID))

  #_______________________________________________________________________________
  # PLOT GENEBODY COVERAGE FOR DEFINE GENE AND SAVE PLOT
  g = ggplot(geneBodyCov, aes(x = Percentiles, y = PercVal*100, color = Condition, linetype = Shape, fill = Condition)) +
    theme_classic() +
    geom_smooth( alpha = .05,se = F) +
    ggtitle(geneOfInterest) +
    scale_color_manual("Samples", values = colors_input) +
    scale_fill_manual("Samples", values = colors_input) +
    scale_x_continuous(breaks = c(0, 50, 100), expand = c(0.01,0.01), labels = c("5'", "mid gene", "3'")) + # change position labels from 0 to 100 to 5' to 3'
    ylab("Relative Coverage (%)") +
    xlab("Gene position") +
    theme(legend.position = "None", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10, color = "black"), 
    axis.text = element_text(color = "black"), axis.ticks.x = element_blank(), 
    panel.grid.major.y = element_line(color = "gray"),
          panel.border = element_rect(colour = "black", fill=NA, size = 1)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,100))

  ggsave(paste0("genebody_Coverage_", geneOfInterest, ".svg"), g_new, width = 10, height = 8, unit = "cm", dpi = 300)
  ggsave(paste0("genebody_Coverage_", geneOfInterest, ".pdf"), g, width = 10, height = 8, unit = "cm", dpi = 300)
}
