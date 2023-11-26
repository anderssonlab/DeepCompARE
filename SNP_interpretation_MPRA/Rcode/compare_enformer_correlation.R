#------------------------------------------------------------------
# Compare various ways to predict effect size using Enformer output
#------------------------------------------------------------------


setwd("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_MPRA")
source("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_MPRA/Rcode/function_utils.R")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/variant_interpretation.R")  
    # aggregate_enformer_delta
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

TRACK_MEANING_FILE <- "/maps/projects/ralab/people/pcr980/Resource_Enformer/targets_human.txt"


#---------------------------------------------------------
# Step 2: calculate prediction and correlation by various methods
#---------------------------------------------------------
# calculate delta of forward and reverse strand
delta_forward <- calculate_enformer_delta("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv",
                                          "Pd3_Enformer_predictions/ref_seq_predictions_forward.csv",
                                          "Pd3_Enformer_predictions/alt_seq_predictions_forward.csv",
                                          "Pd1_MPRA_with_seq_info/ref_seq_forward.csv")
delta_reverse <- calculate_enformer_delta("Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv",
                                          "Pd3_Enformer_predictions/ref_seq_predictions_reverse.csv",
                                          "Pd3_Enformer_predictions/alt_seq_predictions_reverse.csv",
                                          "Pd1_MPRA_with_seq_info/ref_seq_reverse.csv")

# aggregate 5313 predictions and append to df
df <- read.csv("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv")
forward_predictions <- aggregate_enformer_delta(delta_forward,"_forward")
reverse_predictions <- aggregate_enformer_delta(delta_reverse,"_reverse")
fr_predictions <- aggregate_enformer_delta((delta_forward+delta_reverse)/2,"_fr")
df <- cbind(df,forward_predictions,reverse_predictions,fr_predictions)

# calculate prediction correlation by element
prediction_cols <- c(colnames(forward_predictions),colnames(reverse_predictions),colnames(fr_predictions))
correlations <- df %>% 
  group_by(Element) %>%
  summarise(across(all_of(prediction_cols), ~cor(Value, .x)))


#-----------------------------------------
# Step 3: plot and compare between method
#-----------------------------------------
df_plot <- melt(correlations, id.vars = "Element", variable.name = "Method", value.name = "Correlation")
# remove rows with irrelevant cell type

df_plot <- subset(df_plot, !((Element %in% c('F9', 'LDLR', 'LDLR.2','SORT1',"SORT1-flip","SORT1.2")) & grepl("k562", Method)) & 
                          !((Element %in% c('BCL11A','GP1BA', 'HBB', 'HBG1', 'PKLR-24h', 'PKLR-48h')) & grepl("hepg2", Method)))

df_plot$Method <- gsub("hepg2","cell_type_matched",df_plot$Method)
df_plot$Method <- gsub("k562","cell_type_matched",df_plot$Method)


idx <- which(grepl("forward|fr",df_plot$Method) & grepl("dhs",df_plot$Method))
p <- ggplot(df_plot[idx,], aes(x = Element, y = Correlation, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "FR is better than foward",
       x = "Element",
       y = "Correlation Value",
       fill = "Method")+
  scale_fill_brewer(palette="Paired")
print(p)
ggsave(paste0("Generated_plots/Enformer_forward_reverse_and_avg/fr_vs_forward.pdf"), plot = p, width = 10, height = 6)
  


idx <- which(grepl("reverse|fr",df_plot$Method) & grepl("dhs",df_plot$Method))
p <- ggplot(df_plot[idx,], aes(x = Element, y = Correlation, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "FR is better than reverse",
       x = "Element",
       y = "Correlation Value",
       fill = "Method")+
  scale_fill_brewer(palette="Paired")
print(p)
ggsave(paste0("Generated_plots/Enformer_forward_reverse_and_avg/fr_vs_reverse.pdf"), plot = p, width = 10, height = 6)


