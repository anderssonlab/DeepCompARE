setwd("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_raQTL")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/variant_interpretation.R") 
   # aggregate_enformer_delta
source("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_raQTL/Rcode/function_utils.R") 
   # calculate_correlations
   # calculate_enformer_delta_matched_ref_alt
   # calculate_enformer_correlations_avg_fr
   # get_4_correlations_avg_fr

library(reshape2)
library(ggplot2)

TRACK_MEANING_FILE <- "/maps/projects/ralab/people/pcr980/Resource_Enformer/targets_human.txt"
TRUTH_COLUMN <- "log2FC"



calculate_enformer_correlations_f_or_r <- function(pred_ref,pred_alt,truth){
  delta <- calculate_enformer_delta_matched_ref_alt(pred_ref,pred_alt)
  enformer_predictions <- aggregate_enformer_delta(delta,"") 
  calculate_correlations(truth,TRUTH_COLUMN,enformer_predictions)
}


analyze <- function(cell_type){
  cor_forward <- calculate_enformer_correlations_f_or_r(sprintf("Pd3_Enformer_predictions/ref_seq_predictions_%s_sure_raQTL.csv",cell_type),
                                                        sprintf("Pd3_Enformer_predictions/alt_seq_predictions_%s_sure_raQTL.csv",cell_type),
                                                        sprintf("Raw_data/%s_sure_raQTL.csv",cell_type))
  
  cor_reverse <- calculate_enformer_correlations_f_or_r(sprintf("Pd3_Enformer_predictions/ref_seq_predictions_rc_%s_sure_raQTL.csv",cell_type),
                                                        sprintf("Pd3_Enformer_predictions/alt_seq_predictions_rc_%s_sure_raQTL.csv",cell_type),
                                                        sprintf("Raw_data/%s_sure_raQTL.csv",cell_type))
  
  cor_fr <- calculate_enformer_correlations_avg_fr(sprintf("Pd3_Enformer_predictions/ref_seq_predictions_%s_sure_raQTL.csv",cell_type),
                                                   sprintf("Pd3_Enformer_predictions/alt_seq_predictions_%s_sure_raQTL.csv",cell_type),
                                                   sprintf("Pd3_Enformer_predictions/ref_seq_predictions_rc_%s_sure_raQTL.csv",cell_type),
                                                   sprintf("Pd3_Enformer_predictions/alt_seq_predictions_rc_%s_sure_raQTL.csv",cell_type),
                                                   sprintf("Raw_data/%s_sure_raQTL.csv",cell_type))
  
  df_plot <- data.frame(list(
    forward=cor_forward,
    reverse=cor_reverse,
    avg_fr=cor_fr
  ))
  
  df_plot$Method <- rownames(df_plot)
  rownames(df_plot) <- NULL
  df_plot <- melt(df_plot,id.vars = "Method", variable.name = "Strand", value.name = "Correlation")
  
  
  p <- ggplot(df_plot)+
    geom_bar(aes(x=Method,y=Correlation,fill=Strand),position="dodge", stat="identity")+
    ggtitle(cell_type)+
    theme_minimal()
  print(p)
  
  ggsave(paste0("Generated_plots/Enformer_forward_reverse/",cell_type,".pdf"),plot=p,width = 10, height = 6)
  
  
}

analyze("hepg2")
analyze("k562")

