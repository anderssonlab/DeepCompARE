TRACK_MEANING_FILE <- "/maps/projects/ralab/people/pcr980/Resource_Enformer/targets_human.txt"

source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/variant_interpretation.R") 
    # aggregate_enformer_delta


TRUTH_COLUMN <- "log2FC"

calculate_enformer_delta_matched_ref_alt  <- function(df_ref,df_alt){
  # Input:
  #     df_ref: enformer prediction on reference sequence, should have dimension #seqs x 5313. Can be either data.frame or path.
  #     df_alt: enformer prediction on alternative sequence, should have dimension #seqs x 5313. Can be either data.frame or path.
  # Output:
  #.    delta: df_alt-df_ref after log1p transformation on all CAGE tracks
  if (is.character(df_ref)){
    df_ref <- read.csv(df_ref,header = F)
  }
  if (is.character(df_alt)){
    df_alt <- read.csv(df_alt,header = F)
  }
  track_meaning<-read.csv(TRACK_MEANING_FILE,sep="\t", header=TRUE,row.names = 1)
  cage_idx <- grep("CAGE",track_meaning$description)
  df_ref[,cage_idx] <- log1p(df_ref[,cage_idx])
  df_alt[,cage_idx] <- log1p(df_alt[,cage_idx])
  
  delta <- df_alt-df_ref
}


calculate_enformer_correlations_avg_fr <- function(pred_ref1,pred_alt1,pred_ref2,pred_alt2,truth){
  # Input: all can be either data frames or file paths
  #.    pred_ref1: Enformer prediction on reference sequence
  #.    pred_alt1: Enformer prediction on alternative sequence. Should pair up with pred_ref1
  #     pred_ref2: Enformer prediction on reference sequence (reverse complement)
  #     pref_alt2: Enformer prediction on alternative sequence (reverse complement). Should pair up with pred_ref2
  #     truth: truth file
  delta1 <- calculate_enformer_delta_matched_ref_alt(pred_ref1,pred_alt1)
  delta2 <- calculate_enformer_delta_matched_ref_alt(pred_ref2,pred_alt2)
  delta <- (delta1+delta2)/2
  enformer_predictions <- aggregate_enformer_delta(delta,"") 
  calculate_correlations(truth,TRUTH_COLUMN,enformer_predictions)
}


calculate_correlations <- function(df_truth,truth_col,delta){
  # Input:
  #     df_truth: data frame containing truth, or the path
  #.    truth_col: column of truth
  #     delta 
  # Output:
  #    correlation between truth_col and every column of delta
  if (is.character(df_truth)){
    df_truth <- read.csv(df_truth)
  }
  df <- cbind(df_truth,delta)
  df <- df[!is.infinite(df$log2FC),]
  sapply(colnames(delta), function(this_column){
    cor(df[,truth_col],df[,this_column])
  })
}


get_4_correlations_avg_fr <- function(pred_file1,pred_file2,truth_file,truth_col){
  cell_type <- ifelse(grepl("hepg2", truth_file), "hepg2", "k562")
  message("Analyzing ",cell_type)
  delta1 <- calculate_delta(pred_file1,cell_type)
  delta2 <- calculate_delta(pred_file2,cell_type)
  delta <- (delta1+delta2)/2
  colnames(delta) <- c("cage","dhs","starr","sure")
  df_truth <- read.csv(truth_file)
  calculate_correlations(df_truth,TRUTH_COLUMN,delta)
}

