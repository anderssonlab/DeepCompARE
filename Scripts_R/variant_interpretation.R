# This file contains functions to calculate delta for variant interpretation
# For both DeepCompare and Enformer


library(Biostrings)
ENFORMER_TRACK_INFO <- "/maps/projects/ralab/people/pcr980/Resource/Enformer_info/targets_human.txt"



# Input: 
#    required: data frame (or its path) of 16 predictions,  should be direct output of predict_SNP_effect.py
#    optional: cell type in small letter 
# Output: delta 

calculate_delta <- function(df_pred,cell_type=NULL){
  if(is.character(df_pred)){
    df_pred <- read.csv(df_pred)
  }
  
  if(!is.null(cell_type)){
    idx <- grep(cell_type,colnames(df_pred))
    df_pred <- df_pred[,idx]
    delta <- df_pred[,5:8]-df_pred[,1:4]
    colnames(delta) <- c("delta_cage","delta_dhs","delta_starr","delta_sure")
    return(delta)
  }
  delta <- df_pred[,9:16]-df_pred[,1:8]
  colnames(delta) <- gsub("_alt","",colnames(delta))
  colnames(delta) <- gsub("^","delta_",colnames(delta))
  return(delta)
}


# given sequences in character format, write reverse complement in character format
get_rc_from_char_vec <- function(seqs){
  seqs <- DNAStringSet(seqs)
  rev_complements <- reverseComplement(seqs)
  as.character(rev_complements)
}



aggregate_enformer_tracks <- function(enformer_res,suffix=""){
  track_meaning<-read.csv(ENFORMER_TRACK_INFO,sep="\t", header=TRUE,row.names = 1)
  avg <- rowMeans(enformer_res)
  cage_avg <- rowMeans(enformer_res[,grep("CAGE",track_meaning$description)])
  dhs_avg <- rowMeans(enformer_res[,grep("DNASE",track_meaning$description)])
  cage_dhs_avg <- rowMeans(enformer_res[,grep("CAGE|DNASE",track_meaning$description)])
  cage_hepg2_avg <-enformer_res[,which(grepl("HepG2",track_meaning$description) & grepl("CAGE",track_meaning$description))]
  cage_k562_avg <- rowMeans(enformer_res[,grepl("K562",track_meaning$description) & grepl("CAGE",track_meaning$description)])
  dhs_hepg2_avg <- rowMeans(enformer_res[,grepl("HepG2",track_meaning$description) & grepl("DNASE",track_meaning$description)])
  dhs_k562_avg <- rowMeans(enformer_res[,grepl("K562",track_meaning$description) & grepl("DNASE",track_meaning$description)])
  df_summary <- data.frame(cbind(avg,cage_avg,dhs_avg,cage_dhs_avg,cage_hepg2_avg,cage_k562_avg,dhs_hepg2_avg,dhs_k562_avg))
  colnames(df_summary) <- gsub("$",suffix,colnames(df_summary))
  df_summary
}