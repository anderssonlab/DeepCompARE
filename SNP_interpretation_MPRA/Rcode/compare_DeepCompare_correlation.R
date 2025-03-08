#----------------------------------------------------------------------------------------
# Correlate DeepCompare results to experimental value in raw data.
# Compare various ways to calculate predictions (forward, reverse, average)
#---------------------------------------------------------------------------------------


#setwd("/binf-isilon/alab/people/pcr980/SNP_interpretation_MPRA/")
setwd("~/Documents/DeepCompare/SNP_interpretation_MPRA/")
source("~/Documents/DeepCompare/Rcode/variant_interpretation.R") # get calculate_delta()
source("~/Documents/DeepCompare/SNP_interpretation_MPRA/Rcode/function_utils.R")



library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
#-----------------------------------------------
# Step1: use prediction from either forward or reverse
#-----------------------------------------------
correlation_by_element <- function(file_truth,file_pred){
  df_truth <- read.csv(file_truth)
  df_pred <- read.csv(file_pred)
  delta <- calculate_delta(df_pred)
  df <- cbind(df_truth,delta)
  df$X <- NULL
  
  correlations <- df %>% 
    group_by(Element) %>%
    summarise(across(all_of(colnames(delta)), ~cor(Value, .x)))
}

cor_class_forward <- correlation_by_element("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv","Pd2_DeepCompare_predictions/z_forward.csv")
cor_class_reverse <- correlation_by_element("Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv","Pd2_DeepCompare_predictions/z_reverse.csv")
cor_reg_forward <- correlation_by_element("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv","Pd2_DeepCompare_predictions/log1p_signal_forward.csv")
cor_reg_reverse <- correlation_by_element("Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv","Pd2_DeepCompare_predictions/log1p_signal_reverse.csv")

#-----------------------------------------------
# Step 2: use prediction from both forward and reverse
#-----------------------------------------------

cor_class <- correlation_avg_fr_by_element("Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv","Pd2_DeepCompare_predictions/z_forward.csv","Pd2_DeepCompare_predictions/z_reverse.csv")
cor_reg <- correlation_avg_fr_by_element("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv","Pd2_DeepCompare_predictions/log1p_signal_forward.csv","Pd2_DeepCompare_predictions/log1p_signal_reverse.csv")


#-----------------------------------------------
# Step 3: compare the 6 correlations
#-----------------------------------------------
cor_class_long <- melt(cor_class, id.vars = "Element")
cor_class_forward_long <- melt(cor_class_forward, id.vars = "Element")
cor_class_reverse_long <- melt(cor_class_reverse, id.vars = "Element")
cor_reg_long <- melt(cor_reg, id.vars = "Element")
cor_reg_forward_long <- melt(cor_reg_forward, id.vars = "Element")
cor_reg_reverse_long <- melt(cor_reg_reverse, id.vars = "Element")

df_plot <- bind_rows(
  cor_class_long %>% mutate(Dataset = "cor_class"),
  cor_class_forward_long %>% mutate(Dataset = "cor_class_forward"),
  cor_class_reverse_long %>% mutate(Dataset = "cor_class_reverse"),
  cor_reg_long %>% mutate(Dataset = "cor_reg"),
  cor_reg_forward_long %>% mutate(Dataset = "cor_reg_forward"),
  cor_reg_reverse_long %>% mutate(Dataset = "cor_reg_reverse")
  ) %>%
  filter(
    (Element == "F9" & grepl("hepg2", variable)) |
      (Element == "GP1BA" & grepl("k562", variable)) |
      (Element == "HBB" & grepl("k562", variable)) |
      (Element == "HBG1" & grepl("k562", variable)) |
      (grepl("LDLR", Element) & grepl("hepg2", variable)) |
      (grepl("PKLR", Element) & grepl("k562", variable)) |
      (grepl("BCL11A", Element) & grepl("k562", variable)) |
      (grepl("SORT1", Element) & grepl("hepg2", variable))
  )

df_plot$variable <- gsub("_hepg2","",df_plot$variable)
df_plot$variable <- gsub("_k562","",df_plot$variable)
df_plot$variable <- gsub("delta_z_","",df_plot$variable)
df_plot$variable <- gsub("delta_log1p_signal_","",df_plot$variable)

sapply(c("cage","dhs","starr","sure"),function(assay){
  p <- ggplot(df_plot[df_plot$variable==assay,], aes(x=Element, y=value, fill=Dataset)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    theme_minimal() +
    labs(title=assay, x="Element", y="Pearson correlation")+
    scale_fill_brewer(palette="RdBu")
  ggsave(paste0("Generated_plots/DeepCompare_forward_reverse_and_avg/",assay,".pdf"), plot = p, width = 10, height = 6)
})

#--------------------------------------------------------
# Step 4: which assay is most sensitive to strandedness
#--------------------------------------------------------
data <- df_plot %>%
  filter(Dataset %in% c("cor_class_forward","cor_class_reverse","cor_reg_forward","cor_reg_reverse")) %>%
  mutate(type = ifelse(grepl("class", Dataset), "cor_class", "cor_reg"),
         direction = ifelse(grepl("forward", Dataset), "forward", "reverse"))
data$Dataset <- NULL

result <- data %>%
  group_by(Element, variable) %>%
  spread(direction, value) %>%
  mutate(difference = forward - reverse)

p <- ggplot(data=result[result$type=="cor_class",], aes(x=Element,y=difference,fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  labs(title="Classification branch: forward-reverse", y="Correlation difference")+
  theme_minimal()
ggsave(paste0("Generated_plots/DeepCompare_forward_reverse_and_avg/diff_classification.pdf"),plot=p,width = 10, height = 6)

p <- ggplot(data=result[result$type=="cor_reg",], aes(x=Element,y=difference,fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  labs(title="regression branch: forward-reverse", y="Correlation difference")+
  theme_minimal()
ggsave(paste0("Generated_plots/DeepCompare_forward_reverse_and_avg/diff_regression.pdf"),plot=p,width = 10, height = 6)
