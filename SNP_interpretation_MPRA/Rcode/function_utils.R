source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/variant_interpretation.R") 
    # calculate_delta

TRACK_MEANING_FILE <- "/maps/projects/ralab/people/pcr980/Resource_Enformer/targets_human.txt"





correlation_avg_fr_by_element <- function(file_truth,file_pred1,file_pred2){
  # for DeepCompare only
  df_truth <- read.csv(file_truth)
  df_pred1 <- read.csv(file_pred1)
  delta1 <- calculate_delta(df_pred1)
  df_pred2 <- read.csv(file_pred2)
  delta2 <- calculate_delta(df_pred2)
  delta <- (delta1+delta2)/2
  df <- cbind(df_truth,delta)
  df$X <- NULL
  
  correlations <- df %>% 
    group_by(Element) %>%
    summarise(across(all_of(colnames(delta)), ~cor(Value, .x)))
}

calculate_enformer_delta <- function(truth_file, pred_ref, pred_alt,ref_seq_info){
  # read in data
  df <- read.csv(truth_file)
  df_alt <- read.csv(pred_alt,header = F)
  df_ref <- read.csv(pred_ref,header = F)
  
  # log1p transform all CAGE
  track_meaning<-read.csv(TRACK_MEANING_FILE,sep="\t", header=TRUE,row.names = 1)
  cage_idx <- grep("CAGE",track_meaning$description)
  df_ref[,cage_idx] <- log1p(df_ref[,cage_idx])
  df_alt[,cage_idx] <- log1p(df_alt[,cage_idx])
  
  # align df_ref with df_alt
  df_ref_seq_info <- read.csv(ref_seq_info,row.names = 1)
  
  rownames(df_ref) <- df_ref_seq_info$RE_name
  rm(df_ref_seq_info)
  df_ref <- df_ref[df$RE_name,]
  
  # return delta
  df_alt-df_ref
}




