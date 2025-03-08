#-------------------------------------------
# plot DeepCompare and Enformer results.
#-------------------------------------------

library(reshape2)
library(dplyr)
library(ggplot2)
setwd("~/Documents/DeepCompare/SNP_interpretation_MPRA")
source("~/Documents/DeepCompare/SNP_interpretation_MPRA/Rcode/function_utils.R")



# calculate deep compare correlations
cor_dc_class <- correlation_avg_fr_by_element("Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv","Pd2_DeepCompare_predictions/z_forward.csv","Pd2_DeepCompare_predictions/z_reverse.csv")
cor_dc_reg <- correlation_avg_fr_by_element("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv","Pd2_DeepCompare_predictions/log1p_signal_forward.csv","Pd2_DeepCompare_predictions/log1p_signal_reverse.csv")

# calculate enformer correlations
delta_forward <- calculate_enformer_delta("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv",
                                          "Pd3_Enformer_predictions/ref_seq_predictions_forward.csv",
                                          "Pd3_Enformer_predictions/alt_seq_predictions_forward.csv",
                                          "Pd1_MPRA_with_seq_info/ref_seq_forward.csv")
delta_reverse <- calculate_enformer_delta("Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv",
                                          "Pd3_Enformer_predictions/ref_seq_predictions_reverse.csv",
                                          "Pd3_Enformer_predictions/alt_seq_predictions_reverse.csv",
                                          "Pd1_MPRA_with_seq_info/ref_seq_reverse.csv")
enformer_predictions <- aggregate_enformer_delta((delta_forward+delta_reverse)/2,"")
df <- read.csv("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv")
df <- cbind(df,enformer_predictions)
cor_enformer <- df %>% 
  group_by(Element) %>%
  summarise(across(all_of(colnames(enformer_predictions)), ~cor(Value, .x)))



#--------------------------
# Merge correlation values
#--------------------------
cor_dc_class <- melt(cor_dc_class, id.vars = "Element")
colnames(cor_dc_class) <- c("Element","Method","Correlation")
cor_dc_class <- subset(cor_dc_class, !((Element %in% c('F9', 'LDLR', 'LDLR.2','SORT1',"SORT1-flip","SORT1.2")) & grepl("k562", Method)) & 
                    !((Element %in% c('BCL11A','GP1BA', 'HBB', 'HBG1', 'PKLR-24h', 'PKLR-48h')) & grepl("hepg2", Method)))
cor_dc_class$Model <- "DeepCompare_classification"



cor_dc_reg <- melt(cor_dc_reg, id.vars = "Element")
colnames(cor_dc_reg) <- c("Element","Method","Correlation")
cor_dc_reg <- subset(cor_dc_reg, !((Element %in% c('F9', 'LDLR', 'LDLR.2','SORT1',"SORT1-flip","SORT1.2")) & grepl("k562", Method)) & 
                         !((Element %in% c('BCL11A','GP1BA', 'HBB', 'HBG1', 'PKLR-24h', 'PKLR-48h')) & grepl("hepg2", Method)))
# dicide to add "regression" suffix or not
cor_dc_reg$Model <- "DeepCompare_"



cor_enformer <- melt(cor_enformer, id.vars = "Element")
colnames(cor_enformer) <- c("Element","Method","Correlation")
cor_enformer <- subset(cor_enformer, !((Element %in% c('F9', 'LDLR', 'LDLR.2','SORT1',"SORT1-flip","SORT1.2")) & grepl("k562", Method)) & 
                       !((Element %in% c('BCL11A','GP1BA', 'HBB', 'HBG1', 'PKLR-24h', 'PKLR-48h')) & grepl("hepg2", Method)))
cor_enformer <- cor_enformer[cor_enformer$Method %in% c("cage_avg","dhs_avg","cage_hepg2_avg",
                                                        "cage_k562_avg","dhs_hepg2_avg","dhs_k562_avg"),]
cor_enformer$Model <- "Enformer"

# decide to include cor_dc_class or not
correlations <- rbind(cor_dc_reg,cor_enformer)
rownames(correlations) <- NULL

#--------------------------
# Plot correlations
#--------------------------
correlations$Method <- gsub("delta_z_","",correlations$Method)
correlations$Method <- gsub("delta_log1p_signal_","",correlations$Method)
correlations$Method <- gsub("cage_avg","cage_cell_type_agnostic",correlations$Method)
correlations$Method <- gsub("dhs_avg","dhs_cell_type_agnostic",correlations$Method)
correlations$Method <- gsub("hepg2_avg","cell_type_matched",correlations$Method)
correlations$Method <- gsub("k562_avg","cell_type_matched",correlations$Method)
correlations$Method <- gsub("dhs_avg","dhs_cell_type_agnostic",correlations$Method)
correlations$Method <- gsub("_hepg2","",correlations$Method)
correlations$Method <- gsub("_k562","",correlations$Method)
correlations$Method <- gsub("cage$","cage_cell_type_matched",correlations$Method)
correlations$Method <- gsub("dhs$","dhs_cell_type_matched",correlations$Method)
correlations$Method <- gsub("starr$","starr_cell_type_matched",correlations$Method)
correlations$Method <- gsub("sure$","sure_cell_type_matched",correlations$Method)
correlations$Model <- gsub("DeepCompare_", "DeepCompare", correlations$Model)



correlations <- correlations[!correlations$Element%in% c("LDLR.2","PKLR-48h","SORT1-flip","SORT1.2"),]
correlations$Element <- gsub("PKLR-24h","PKLR",correlations$Element)


# Calculate medians
medians <- correlations %>%
  group_by(Element, Model) %>%
  summarize(median_correlation = median(Correlation, na.rm = TRUE))

# Merge the medians back into the original data
correlations_with_medians <- left_join(correlations, medians, by = c("Element", "Model"))


p <- ggplot(correlations_with_medians, aes(x = Model, y = Correlation, color=Model, shape=Method)) + 
  geom_errorbar(aes(ymin = median_correlation, ymax = median_correlation), 
                width = 0.2, position = position_dodge(width = 0.2),linetype = "solid",color="black") +
  geom_point(position = position_dodge(width = 0.5),size=3) + 
  scale_shape_manual(values = c(15,16,17,18,4,8))+
  facet_wrap(~ Element, scales = "free_y", ncol = 2) +  
  theme(axis.text.x = element_blank(), 
        axis.text.x.bottom = element_text(),
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove panel background
        plot.background = element_blank(),
        #panel.grid.major = element_line(color = "grey"),
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 10)
        ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)), # Increase size of color legend keys
    shape = guide_legend(override.aes = list(size = 3))
  )+
  labs(title = "MPRA Variant interpretation: DeepCompare v.s. Enformer", x = "Model", y = "Pearson Correlation")

ggsave(paste0("Generated_plots/Publication/MPRA_variant_interpretation_DeepCompare_vs_Enformer.pdf"), plot = p, width = 10, height = 10)


