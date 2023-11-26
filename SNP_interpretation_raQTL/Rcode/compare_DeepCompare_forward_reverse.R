setwd("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_raQTL")
source("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_raQTL/Rcode/function_utils.R") 
    # calculate_correlations
    # get_4_correlations_avg_fr
library(reshape2)
library(ggplot2)




get_4_correlations_f_or_r <- function(pred_file,truth_file,truth_col){
  cell_type <- ifelse(grepl("hepg2", truth_file), "hepg2", "k562")
  message("Analyzing ",cell_type)
  delta <- calculate_delta(pred_file,cell_type)
  colnames(delta) <- c("cage","dhs","starr","sure")
  
  df_truth <- read.csv(truth_file)
  calculate_correlations(df_truth,truth_col,delta)
}



infer_file_names <- function(cell_type,model){
  prefix <- ifelse(model=="classification","z","log1p_signal")
  pred_file1 <- sprintf("Pd2_DeepCompare_predictions/%s_%s_forward.csv",prefix,cell_type)
  pred_file2 <- sprintf("Pd2_DeepCompare_predictions/%s_%s_reverse.csv",prefix,cell_type)
  truth_file <- sprintf("Raw_data/%s_sure_raQTL.csv",cell_type)
  return(c(pred_file1,pred_file2,truth_file))
}

analyze <- function(cell_type,model,truth_col="log2FC"){
  inferred_names <- infer_file_names(cell_type,model)
  pred_file1 <- inferred_names[1]
  pred_file2 <- inferred_names[2]
  truth_file <- inferred_names[3]
  message("Processing ", pred_file1, "\n", pred_file2, "\n",truth_file )
  cor_forward <- get_4_correlations_f_or_r(pred_file1,truth_file,truth_col)
  cor_reverse <- get_4_correlations_f_or_r(pred_file2,truth_file,truth_col)
  cor_fr <- get_4_correlations_avg_fr(pred_file1,pred_file2,truth_file,truth_col)
  
  df_plot <- data.frame(list(
    forward=cor_forward,
    reverse=cor_reverse,
    avg_fr=cor_fr
  ))
  df_plot$assay <- c("cage","dhs","starr","sure")
  df_plot <- melt(df_plot,id.vars = "assay", variable.name = "Method", value.name = "Correlation")
  
  
  p <- ggplot(df_plot)+
    geom_bar(aes(x=assay,y=Correlation,fill=Method),position="dodge", stat="identity")+
    ggtitle(paste(cell_type,model))+
    theme_minimal()
  print(p)
  
  ggsave(paste0("Generated_plots/DeepCompare_forward_reverse/",cell_type,"_",model,".pdf"),plot=p,width = 10, height = 6)
  
}

analyze("hepg2","classification")
analyze("hepg2","regression")
analyze("k562","classification")
analyze("k562","regression")

