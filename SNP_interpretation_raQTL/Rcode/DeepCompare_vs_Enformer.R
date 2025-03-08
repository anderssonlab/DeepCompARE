setwd("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_raQTL")
source("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_raQTL/Rcode/function_utils.R")
    # calculate_correlations
    # calculate_enformer_correlations_avg_fr
    # get_4_correlations_avg_fr

library(ggplot2)
library(dplyr)


cor_enformer_hepg2 <- calculate_enformer_correlations_avg_fr("Pd3_Enformer_predictions/ref_seq_predictions_hepg2_sure_raQTL.csv",
                                                             "Pd3_Enformer_predictions/alt_seq_predictions_hepg2_sure_raQTL.csv",
                                                             "Pd3_Enformer_predictions/ref_seq_predictions_rc_hepg2_sure_raQTL.csv",
                                                             "Pd3_Enformer_predictions/alt_seq_predictions_rc_hepg2_sure_raQTL.csv",
                                                             "Raw_data/hepg2_sure_raQTL.csv")

cor_enformer_k562 <- calculate_enformer_correlations_avg_fr("Pd3_Enformer_predictions/ref_seq_predictions_k562_sure_raQTL.csv",
                                                            "Pd3_Enformer_predictions/alt_seq_predictions_k562_sure_raQTL.csv",
                                                            "Pd3_Enformer_predictions/ref_seq_predictions_rc_k562_sure_raQTL.csv",
                                                            "Pd3_Enformer_predictions/alt_seq_predictions_rc_k562_sure_raQTL.csv",
                                                            "Raw_data/k562_sure_raQTL.csv")

cor_deepcompare_class_hepg2 <- get_4_correlations_avg_fr("Pd2_DeepCompare_predictions/z_hepg2_forward.csv",
                                                         "Pd2_DeepCompare_predictions/z_hepg2_reverse.csv",
                                                         "Raw_data/hepg2_sure_raQTL.csv",
                                                         "log2FC")

cor_deepcompare_reg_hepg2 <- get_4_correlations_avg_fr("Pd2_DeepCompare_predictions/log1p_signal_hepg2_forward.csv",
                                                       "Pd2_DeepCompare_predictions/log1p_signal_hepg2_reverse.csv",
                                                       "Raw_data/hepg2_sure_raQTL.csv",
                                                       "log2FC")

cor_deepcompare_class_k562 <- get_4_correlations_avg_fr("Pd2_DeepCompare_predictions/z_k562_forward.csv",
                                                         "Pd2_DeepCompare_predictions/z_k562_reverse.csv",
                                                         "Raw_data/k562_sure_raQTL.csv",
                                                         "log2FC")

cor_deepcompare_reg_k562 <- get_4_correlations_avg_fr("Pd2_DeepCompare_predictions/log1p_signal_k562_forward.csv",
                                                       "Pd2_DeepCompare_predictions/log1p_signal_k562_reverse.csv",
                                                       "Raw_data/k562_sure_raQTL.csv",
                                                       "log2FC")


cor_enformer_hepg2 <- cor_enformer_hepg2[c(2,3,5,7)]
cor_enformer_k562 <- cor_enformer_k562[c(2,3,6,8)]

methods <- c(names(cor_enformer_hepg2),
             names(cor_enformer_k562),
             names(cor_deepcompare_class_hepg2),
             names(cor_deepcompare_reg_hepg2),
             names(cor_deepcompare_class_k562),
             names(cor_deepcompare_reg_k562))
correlations <- c(unname(cor_enformer_hepg2),
                  unname(cor_enformer_k562),
                  unname(cor_deepcompare_class_hepg2),
                  unname(cor_deepcompare_reg_hepg2),
                  unname(cor_deepcompare_class_k562),
                  unname(cor_deepcompare_reg_k562))
models <-c(rep("Enformer",8),
           rep("DeepCompare\nclassification",4),
           rep("DeepCompare\nregression",4),
           rep("DeepCompare\nclassification",4),
           rep("DeepCompare\nregression",4))

cell_types <- c(rep("HepG2",4),rep("K562",4),
               rep("HepG2",8),rep("K562",8))

df_plot <- data.frame(list(
  Method=methods,
  Correlation=correlations,
  Model=models,
  Cell_type=cell_types
))

df_plot$Method <- gsub("^cage$","cage_cell_type_matched",df_plot$Method)
df_plot$Method <- gsub("^dhs$","dhs_cell_type_matched",df_plot$Method)
df_plot$Method <- gsub("^starr$","starr_cell_type_matched",df_plot$Method)
df_plot$Method <- gsub("^sure$","sure_cell_type_matched",df_plot$Method)
df_plot$Method <- gsub("^cage_avg$","cage_cell_type_agnostic",df_plot$Method)
df_plot$Method <- gsub("^dhs_avg$","dhs_cell_type_agnostic",df_plot$Method)
df_plot$Method <- gsub("hepg2_avg$","cell_type_matched",df_plot$Method)
df_plot$Method <- gsub("k562_avg$","cell_type_matched",df_plot$Method)
df_plot$Cell_type <- gsub("HepG2", "HepG2 raQTL",df_plot$Cell_type)
df_plot$Cell_type <- gsub("K562", "K562 raQTL",df_plot$Cell_type)



medians <- df_plot %>%
  group_by(Cell_type, Model) %>%
  summarize(median_correlation = median(Correlation, na.rm = TRUE))

# Merge the medians back into the original data
df_plot_with_medians <- left_join(df_plot, medians, by = c("Cell_type", "Model"))



p <- ggplot(df_plot_with_medians,(aes(x=Model,y=Correlation,color=Model,shape=Method)))+
  geom_point(position = position_dodge(width = 0.5),size=3)+
  geom_errorbar(aes(ymin = median_correlation, ymax = median_correlation), 
                width = 0.2, position = position_dodge(width = 0.2),linetype = "solid",color="black") +
  scale_shape_manual(values = c(15,16,17,18,4,8))+
  facet_wrap(~ Cell_type, scales = "free_y", ncol = 2)+
  theme(axis.text.x = element_blank(), 
        axis.text.x.bottom = element_text(),
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove panel background
        plot.background = element_blank(),
        #panel.grid.major = element_line(color = "grey"),
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 10)
  ) +
  labs(title = "SuRE-Seq raQTL variant interpretation: DeepCompare v.s. Enformer", 
       x = "Model", 
       y = "Pearson Correlation")

ggsave(paste0("Generated_plots/Publication/raQTL_variant_interpretation_DeepCompare_vs_Enformer.pdf"), plot = p, width = 10, height = 4)


