setwd("~/Documents/DeepCompare/SNP_interpretation_MPRA")
source("~/Documents/DeepCompare/Rcode/plot_importance_by_position.R")
source("~/Documents/DeepCompare/Rcode/variant_interpretation.R")
library(dplyr)
library(tidyr)

prepare_df_plot <- function(df_truth,delta,cell_type){
  delta <- delta[,grepl(cell_type,colnames(delta))]
  df <- cbind(df_truth,delta)
  # reshape data frame
  df_long <- pivot_longer(df, cols = c(Value, colnames(delta)), 
                          names_to = "Measurement", values_to = "Importance")
  
  # rename
  df_long$Measurement <- gsub("Value","Experimental effect size",df_long$Measurement)
  df_long$Measurement <- gsub("delta_z_","",df_long$Measurement)
  df_long$Measurement <- gsub("delta_log1p_signal_","",df_long$Measurement)
  df_long$Measurement <- gsub("_hepg2","",df_long$Measurement)
  df_long$Measurement <- gsub("_k562","",df_long$Measurement)
  df_long$Measurement <- gsub("cage","CAGE prediction",df_long$Measurement)
  df_long$Measurement <- gsub("dhs","DHS prediction",df_long$Measurement)
  df_long$Measurement <- gsub("starr","STARR prediction",df_long$Measurement)
  df_long$Measurement <- gsub("sure","SuRE prediction",df_long$Measurement)
  
  # reorder
  df_long$Measurement <- factor(df_long$Measurement, levels = c("Experimental effect size", "CAGE prediction","DHS prediction", "STARR prediction", "SuRE prediction"))
  df_long
}


plot_importance_by_position_multi_tracks <- function(df_truth,delta,cell_type,title){
  df_plot <- prepare_df_plot(df_truth,delta,cell_type)

  p <- ggplot(df_plot,aes(x=relevant_position, y=Importance)) +
    geom_point(aes(x=relevant_position, y=Importance,color=Alt),size=0.3)+
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_bar(stat="identity",position=position_identity(),width = 0.8) +
    facet_wrap(~Measurement,ncol = 1,scales="free_y") +
    labs(title = title)+
    xlab("Relevant Position") +
    ylab("Effect Size") +
    guides(color = guide_legend(override.aes = list(size=2)))+
    labs(color='Alternative base') +
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()
    )
  print(p)
  ggsave(paste0("Generated_plots/Publication/",title,".pdf"), plot = p, width = 8, height = 8)
}



plot_importance_by_position_multi_tracks_with_boxes <- function(df_truth,delta,cell_type,box_list,title){
  df_plot <- prepare_df_plot(df_truth,delta,cell_type)
  
  p <- ggplot(df_plot,aes(x=relevant_position, y=Importance)) +
    geom_point(aes(x=relevant_position, y=Importance,color=Alt),size=0.3)+
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_bar(stat="identity",position=position_identity(),width = 0.8) +
    facet_wrap(~Measurement,ncol = 1,scales="free_y") +
    labs(title = title)+
    xlab("Relevant Position") +
    ylab("Effect Size") +
    guides(color = guide_legend(override.aes = list(size=2)))+
    labs(color='Alternative base') +
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()
    )
  for (this_box in box_list){
    p <- p+geom_rect(xmin = this_box[1], xmax = this_box[2], ymin = -Inf, ymax = Inf,fill = "cyan3",alpha = 0.01)
  }
  print(p)
  ggsave(paste0("Generated_plots/Publication/",title,".pdf"), plot = p, width = 8, height = 8)
  
  
}








# experimental value
df_truth <- read.csv("Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv")
df_truth$X <- NULL

# prediction, averaging both forward and reverse
df_pred_forward <- read.csv("Pd2_DeepCompare_predictions/z_forward.csv")
delta_forward <- calculate_delta(df_pred_forward)
df_pred_reverse <- read.csv("Pd2_DeepCompare_predictions/z_reverse.csv")
delta_reverse <- calculate_delta(df_pred_reverse)
delta <- (delta_forward+delta_reverse)/2


# plot
idx <- which(df_truth$Element=="F9")
plot_importance_by_position_multi_tracks_with_boxes(df_truth[idx,],delta[idx,],"hepg2",
                                                    list(c(428,435),c(414,419),c(396,404),c(387,401)),
                                                    "Effect size prediction on F9 promoter"
                                                    )

idx <- which(df_truth$Element=="LDLR")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"hepg2","Effect size prediction on LDLR promoter")

idx <- which(df_truth$Element=="PKLR-24h")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"k562","Effect size prediction on PKLR promoter")

idx <- which(df_truth$Element=="SORT1")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"hepg2","Effect size prediction on SORT1 enhancer")

idx <- which(df_truth$Element=="BCL11A")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"k562","Effect size prediction on BCL11A enhancer")

idx <- which(df_truth$Element=="HBB")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"k562","Effect size prediction on HBB promoter")

idx <- which(df_truth$Element=="HBG1")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"k562","Effect size prediction on HBG1 promoter")

idx <- which(df_truth$Element=="GP1BA")
plot_importance_by_position_multi_tracks(df[idx,],delta[idx,],"k562","Effect size prediction on HBG1 promoter")
