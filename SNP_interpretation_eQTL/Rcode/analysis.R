setwd("~/Documents/DeepCompare/eQTL")
source("~/Documents/DeepCompare/Rcode/variant_interpretation.R")


library(ggplot2)
library(patchwork)
library(gridExtra)

# enrichment given z threshold
calculate_odd_ratio <-function(thresh,pip_vec,delta_vec){
  high_z_high_pip <- sum(pip_vec>0.9 & delta_vec>=thresh)
  high_z_low_pip <- sum(pip_vec<0.01 & delta_vec>=thresh)
  low_z_high_pip <- sum(pip_vec>0.9 & delta_vec<thresh)
  low_z_low_pip <- sum(pip_vec<0.01 & delta_vec<thresh)
  odd_high_z <- high_z_high_pip/high_z_low_pip
  odd_low_z <- low_z_high_pip/low_z_low_pip
  if (high_z_low_pip==0) odd_high_z <-0
  if (low_z_low_pip==0) odd_low_z <-0
  odd_ratio <- odd_high_z/odd_low_z
  if (odd_low_z==0) odd_ratio <- 0
  odd_ratio
} 
#fisher exact test/chi-square test, get p values.

# enrichment for every percentile as threshold
or_along_quantile <- function(df,delta,assay){
  quantiles <- quantile(abs(delta[[assay]]), probs = seq(0, 0.99, by = 0.01))
  sapply(quantiles,function(thresh){
    calculate_odd_ratio(thresh,df$pip,abs(delta[[assay]]))
  })
}

# Plot enrichment for one dataset
plot_enrichment <- function(df,df_pred,cell_type,title){
  delta <- calculate_delta(df_pred,cell_type)
  enrichment_cage <- or_along_quantile(df,delta,"cage")
  enrichment_dhs <- or_along_quantile(df,delta,"dhs")
  enrichment_starr <- or_along_quantile(df,delta,"starr")
  enrichment_sure <- or_along_quantile(df,delta,"sure")
  enrichment <- c(enrichment_cage,enrichment_dhs,enrichment_starr,enrichment_sure)
  quantile <- rep(seq(0,99,by=1),4)
  assay <- rep(c("cage","dhs","starr","sure"),each=100)
  df_plot <- data.frame(cbind(enrichment,quantile,assay))
  df_plot$enrichment <- as.numeric(df_plot$enrichment)
  df_plot$quantile <- as.numeric(df_plot$quantile)
  ggplot(df_plot, aes(x=quantile,y=enrichment,color=assay))+
    geom_line()+
    theme_bw()+
    xlab("abs(delta prediction) quantile")+
    ggtitle(title)
}


# Plot enrichment for both regression and classification values
analyze_dataset_by_enrichment <- function(suffix){
  df <- read.table(paste0("Raw_data/QTD000",suffix,".credible_sets.tsv"),header=T)
  df_regression <- read.csv(paste0("Processed_data/log1p_signal_QTD000",suffix,".csv"))
  df_classification <- read.csv(paste0("Processed_data/z_QTD000",suffix,".csv"))
  idx <- which(df$pvalue<0.05)
  df <- df[idx,]
  df_regression <- df_regression[idx,]
  df_classification <- df_classification[idx,]
  p1 <- plot_enrichment(df,df_regression,"hepg2","regression hepg2")
  p2 <- plot_enrichment(df,df_classification,"hepg2","classification hepg2")
  p3 <- plot_enrichment(df,df_regression,"k562","regression k562")
  p4 <- plot_enrichment(df,df_classification,"k562","classification k562")
  pdf(paste0("Generated_plots/",suffix,".pdf"), width = 8, height = 6)
  grid.arrange(p1, p2, p3, p4,ncol = 2)
  dev.off()
}

# scatter plot of delta-z v.s. PIP
plot_scatter <- function(df,df_pred,cell_type,title){
  delta <- calculate_delta(df_pred,cell_type)
  delta$pip <- df$pip
  p1 <- ggplot(data=delta,aes(x=pip,y=cage))+geom_point()+theme_bw()
  p2 <- ggplot(data=delta,aes(x=pip,y=dhs))+geom_point()+theme_bw()
  p3 <- ggplot(data=delta,aes(x=pip,y=starr))+geom_point()+theme_bw()
  p4 <- ggplot(data=delta,aes(x=pip,y=sure))+geom_point()+theme_bw()
  grid.arrange(p1, p2, p3, p4,ncol = 2,top=title)
}


analyze_dataset_by_scatter <- function(suffix,cell_type){
  df <- read.table(paste0("Raw_data/QTD000",suffix,".credible_sets.tsv"),header=T)
  df_regression <- read.csv(paste0("Processed_data/log1p_signal_QTD000",suffix,".csv"))
  df_classification <- read.csv(paste0("Processed_data/z_QTD000",suffix,".csv"))
  
  # subset to contain only significant
  idx <- which(df$pvalue<0.05)
  df_sub <- df[idx,]
  df_regression_sub <- df_regression[idx,]
  df_classification_sub <- df_classification[idx,]

  # start plotting
  all_grobs <- list(plot_scatter(df,df_regression,cell_type, "regression whole"),
                    plot_scatter(df_sub,df_regression_sub,cell_type, "regression significant"),
                    plot_scatter(df,df_classification,cell_type, "classification whole"),
                    plot_scatter(df_sub,df_classification_sub,cell_type, "classification significant")
  )
  
  pdf(paste0("Generated_plots/",suffix,"_scatter.pdf"), width = 8, height = 8)
  grid.arrange(grobs=all_grobs, ncol=2)
  dev.off()
}

#----------------------------------
# Main analaysis
# 404,406,406 are hepatocyte
# 549, 550, 551 are whole blood
#----------------------------------

# enrichment plots
analyze_dataset_by_enrichment(404)
analyze_dataset_by_enrichment(405)
analyze_dataset_by_enrichment(406)

analyze_dataset_by_enrichment(549)
analyze_dataset_by_enrichment(550)
analyze_dataset_by_enrichment(551)

analyze_dataset_by_scatter(404,"hepg2")
analyze_dataset_by_scatter(405,"hepg2")
analyze_dataset_by_scatter(406,"hepg2")

analyze_dataset_by_scatter(549,"k562")
analyze_dataset_by_scatter(550,"k562")
analyze_dataset_by_scatter(551,"k562")

df <-read.csv("Raw_data/QTD000404.credible_sets.tsv",sep="\t")
hist(df$pip)
sum(df$pip>0.9)
df_pred <- read.csv("Processed_data/log1p_signal_QTD000404.csv")
delta <- calculate_delta(df_pred,"hepg2")

plot(df$pip,df$beta)
plot(abs(delta$cage),df$beta)
