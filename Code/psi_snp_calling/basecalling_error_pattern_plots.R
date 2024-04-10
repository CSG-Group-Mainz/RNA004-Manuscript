################################################################################
####                                                                        ####
####       ANALYSIS OF BASECALLING ERROR PATTERNS OF VECTOR SEQUENCES       ####
####                                                                        ####
################################################################################

#_______________________________________________________________________________
# LIBRARIES ####
library(readr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggtext)
library(ggh4x)


#_______________________________________________________________________________
# PLOTTING SCHEME ####
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"), 
             axis.text = element_text(color = "black"),
             strip.text = element_markdown())
ref_cols <- c("EGFP" = "#7fb069", 
              "mCherry" = "#ca3c25")
color_code = c("RNA002_HEK293" = "#0a9396", 
                "RNA004_HEK293" = "#ca6702")
#_______________________________________________________________________________
# FUNCTIONS ####
processMPILEUP_splint <- function(splint_bam){
  splint_bam = as.data.frame(splint_bam)
  pattern = str_replace_all(splint_bam$PAT, "\\^[ -~]", "S")
  splint_bam$A = str_count(toupper(pattern), "A") + ifelse(splint_bam$REF_BASE == "A", str_count(pattern, "\\."), 0)
  splint_bam$T = str_count(toupper(pattern), "T") + ifelse(splint_bam$REF_BASE == "T", str_count(pattern, "\\."), 0)
  splint_bam$G = str_count(toupper(pattern), "G")+ ifelse(splint_bam$REF_BASE == "G", str_count(pattern, "\\."), 0)
  splint_bam$C = str_count(toupper(pattern), "C")+ ifelse(splint_bam$REF_BASE == "C", str_count(pattern, "\\."), 0)
  splint_bam$MATCH = str_count(pattern, "\\.")
  splint_bam$MISMATCH = str_count(pattern, "[ATGCatgc]")
  splint_bam$DELETION = str_count(pattern, "\\*")
  splint_bam$DELETION_after = str_count(pattern, "\\-")
  splint_bam$INSERTION = str_count(pattern, "\\.\\+")
  
  splint_bam$REVERSE = str_count(pattern, "\\,")
  splint_bam$REF_SKIP = str_count(pattern, ">")
  return(splint_bam)
}

process_mpileup <- function(mpileup_splint, sample_levels){
  A_MPILE = mpileup_splint[, c(1:6)]
  B_MPILE = mpileup_splint[, c(1:3, 7:9)]
  colnames(A_MPILE) = c("REF", "POS", "REF_BASE", "COV", "PAT", "QUAL")
  colnames(B_MPILE) = c("REF", "POS", "REF_BASE", "COV", "PAT", "QUAL")
  A_MPILE = processMPILEUP_splint(splint_bam = A_MPILE)
  B_MPILE = processMPILEUP_splint(splint_bam = B_MPILE)
  A_MPILE$MOD = sample_levels[1]
  B_MPILE$MOD = sample_levels[2]
  MPILE = rbind(A_MPILE, B_MPILE)
  MPILE$MATCH = MPILE$COV-(MPILE$MISMATCH+MPILE$DELETION+MPILE$REVERSE+MPILE$REF_SKIP+MPILE$INSERTION)
  MPILE$SUM = MPILE$MATCH + MPILE$INSERTION+MPILE$MISMATCH+MPILE$DELETION+MPILE$REVERSE+MPILE$REF_SKIP
  MPILE$DIFF = MPILE$COV-MPILE$SUM
  print(table(MPILE$DIFF))
  
  MPILE_stack = MPILE %>% 
    gather(key = "condition", value = "count", MATCH, INSERTION, MISMATCH, DELETION, REVERSE, REF_SKIP)
  MPILE_stack$condition[MPILE_stack$condition == "MISMATCH"] = ifelse(MPILE_stack[MPILE_stack$condition == "MISMATCH", "REF_BASE"] %in% c("A", "G", "T"), "MISMATCH_C", "MISMATCH_other")
  
  MPILE_stack$condition_prob = MPILE_stack$count/MPILE_stack$COV
  MPILE_stack$condition = factor(MPILE_stack$condition, levels = rev(c("MISMATCH_C", "MISMATCH_other", "INSERTION", "DELETION", "REVERSE", "REF_SKIP", "MATCH")))
  
  MPILE_stack$MOD = factor(MPILE_stack$MOD, levels = sample_levels)
  
  
  return(MPILE_stack)
}

process_small_window <- function(MPILE_stack){
  MPILE_stack_sub = MPILE_stack[(MPILE_stack$REF == "EGFP" & MPILE_stack$POS >= 110 & MPILE_stack$POS <= 120) | 
                                  (MPILE_stack$REF == "mCherry" & MPILE_stack$POS >= 560 & MPILE_stack$POS <= 570),]
  
  new_labels = rep(c(paste0("-", 5:1), "0", paste0("+",1:5)),2)
  names(new_labels) = c(110:120, 560:570)
  
  MPILE_stack_sub$REL_POS = as.character(new_labels[as.character(MPILE_stack_sub$POS)])
  MPILE_stack_sub$REL_POS = factor(MPILE_stack_sub$REL_POS, c(paste0("-", 5:1), "0", paste0("+",1:5)))
  MPILE_stack_sub$condition_prob = MPILE_stack_sub$count/MPILE_stack_sub$COV
  
  MPILE_stack_sub$condition = as.character(MPILE_stack_sub$condition)
  MPILE_stack_sub_sub = MPILE_stack_sub[-which(MPILE_stack_sub$condition %in% c("MATCH", "REVERSE", "REF_SKIP")),]
  MPILE_stack_sub_sub$condition = factor(MPILE_stack_sub_sub$condition, levels = rev(c("MISMATCH_C", "MISMATCH_other", "INSERTION", "DELETION")))
  MPILE_stack_sub_sub_sub = unique(MPILE_stack_sub_sub[MPILE_stack_sub_sub$condition == "MISMATCH_C", c("REL_POS", "condition_prob", "condition", "MOD", "REF")])
  
  return(MPILE_stack_sub_sub)
}

process_large_window <- function(MPILE_stack){
  MPILE_stack_sub = MPILE_stack[(MPILE_stack$REF == "EGFP" & MPILE_stack$POS >= (115-49) & MPILE_stack$POS <= (115+49)) | 
                                  (MPILE_stack$REF == "mCherry" & MPILE_stack$POS >= (565-49) & MPILE_stack$POS <= (565+49)),]
  
  MPILE_stack_sub$condition = as.character(MPILE_stack_sub$condition)
  MPILE_stack_sub_sub = MPILE_stack_sub[-which(MPILE_stack_sub$condition %in% c("MATCH", "REVERSE", "REF_SKIP")),]
  MPILE_stack_sub_sub$condition = factor(MPILE_stack_sub_sub$condition, levels = rev(c("MISMATCH_C", "MISMATCH_other", "INSERTION", "DELETION")))
  tmp = MPILE_stack_sub_sub %>% group_by(POS, MOD, REF) %>% summarise(n = sum(condition_prob))
  return(tmp)
}

#________________________________
# LOAD INPUT MPILEUP ####
RNA002_Q12_dorado_mpileup <- fread("RNA002_HEK293.q12.mpileup")
RNA004_Q12_dorado_mpileup <- fread("RNA004_HEK293.q12.mpileup")

# DEFINE COLORS FOR BASECALLING ERRORS
my_BC_cols = c("MISMATCH_C" = "#ee6c4e", "MISMATCH_other" = "#102442", "INSERTION" = "#457b9d", "DELETION" = "#98c1d8")

# EXTRACT COUNTS FOR EACH BASECALLING ERROR FROM MPILEUP FILE
MPILE_STACK_RNA004 = process_mpileup(RNA004_Q12_dorado_mpileup, sample_levels = c("Sample_A", "Sample_B"))
MPILE_STACK_RNA002 = process_mpileup(RNA002_Q12_dorado_mpileup, sample_levels = c("Sample_A", "Sample_B"))

# EXTRACT BASECALLING ERROR FOR +/-10 BP WINDOW AROUND MODIFICATION
MPILE_STACK_RNA004_small = process_small_window(MPILE_stack = MPILE_STACK_RNA004)
MPILE_STACK_RNA002_small = process_small_window(MPILE_STACK_RNA002)
MPILE_STACK_RNA004_small$source = "RNA004_HEK293"
MPILE_STACK_RNA002_small$source = "RNA002_HEK293"

MPILE_STACK_small = rbind(MPILE_STACK_RNA002_small, MPILE_STACK_RNA004_small)

# PLOT 20BP WINDOW
g11 = ggplot(MPILE_STACK_small, aes(x = REL_POS, y = condition_prob*100, fill = condition)) +
  geom_bar(stat = "identity", colour = "black") +
  facet_grid(source~MOD+REF) +
  theme_bw() +
  xlab("Relative vector position") +
  ylab("Basecalling error %") +
  #theme(strip.text = element_markdown()) +
  scale_fill_manual("", values = my_BC_cols) + 
  theme(legend.position = "bottom", axis.ticks.x = element_blank(),panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill="white"))

ggsave("Basecalling_Errors_20bp.png", g11, width = 26, height = 10, units = "cm", dpi = 300)
ggsave("Basecalling_Errors_20bp.svg", g11, width = 26, height = 10, units = "cm", dpi = 300)
ggsave("Basecalling_Errors_20bp.pdf", g11, width = 26, height = 10, units = "cm", dpi = 300)

# EXTRACT BASECALLING ERROR FOR +/-50 BP WINDOW AROUND MODIFICATION
MPILE_STACK_RNA004_large = process_large_window(MPILE_stack = MPILE_STACK_RNA004)
MPILE_STACK_RNA002_large = process_large_window(MPILE_stack = MPILE_STACK_RNA002)
MPILE_STACK_RNA004_large$source = "RNA004_HEK293"
MPILE_STACK_RNA002_large$source = "RNA002_HEK293"

MPILE_STACK_large = rbind(MPILE_STACK_RNA002_large, MPILE_STACK_RNA004_large)
# Define x axis scales
df_scales <- data.frame(
  Panel = c("RNA002_blood_A", "RNA002_blood_B"),
  ymin = c(70, 520),
  ymax = c(160, 610),
  n = c(115, 565)
)
df_scales <- split(df_scales, df_scales$Panel)
scales <- lapply(df_scales, function(x) {
  scale_x_continuous(breaks = c(x$ymin,x$n, x$ymax),  labels = c("<- 5'", "PsU", "3' ->"))
})

# PLOT 100BP WINDOW
g = ggplot(MPILE_STACK_large, aes(x = POS, y = n*100, color = source)) +
  facet_grid(MOD~REF, scales = "free") +
  geom_line() +
  geom_area(position=position_identity(),alpha = 0.7, aes(fill = source)) +
  theme_bw() +
  xlab("Relative vector position") +
  ylab("Basecalling error %") +
  scale_fill_manual(values = color_code[c("RNA002_HEK293", "RNA004_HEK293")]) +
  scale_color_manual(values = color_code[c("RNA002_HEK293", "RNA004_HEK293")]) +
  ggh4x::facetted_pos_scales(
    x = scales) + 
  theme(legend.position = "bottom", axis.ticks.x = element_blank(), 
             panel.grid = element_blank(), strip.background = element_rect(fill="white"))
ggsave("Basecalling_Errors_100bp.png", g, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Basecalling_Errors_100bp.svg", g, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Basecalling_Errors_100bp.pdf", g, width = 20, height = 10, units = "cm", dpi = 300)

