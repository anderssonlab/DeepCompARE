setwd("/maps/projects/ralab/people/pcr980/DeepCompare/MPRA_combinatorics")

source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/variant_interpretation.R")


enformer_pred <- read.csv("Pd3_Enformer_predictions/enformer_predictions.csv",header=FALSE)

enformer_pred_aggregated <- aggregate_enformer_tracks(enformer_pred)
write.csv(enformer_pred_aggregated,"Pd4_Enformer_predictions_aggregated/enformer_predictions_aggregated.csv")
