

TRACK_MEANING_FILE <- "~/Documents/DeepCompare/DCMain_Sup2_benchmarkl_enformer_MPRA/Raw_data/enformer_targets_human.txt"




library(Biostrings)



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




aggregate_enformer_delta <- function(delta,suffix){
  track_meaning<-read.csv(TRACK_MEANING_FILE,sep="\t", header=TRUE,row.names = 1)
  avg <- rowMeans(delta)
  cage_avg <- rowMeans(delta[,grep("CAGE",track_meaning$description)])
  dhs_avg <- rowMeans(delta[,grep("DNASE",track_meaning$description)])
  cage_dhs_avg <- rowMeans(delta[,grep("CAGE|DNASE",track_meaning$description)])
  cage_hepg2_avg <-delta[,which(grepl("HepG2",track_meaning$description) & grepl("CAGE",track_meaning$description))]
  cage_k562_avg <- rowMeans(delta[,grepl("K562",track_meaning$description) & grepl("CAGE",track_meaning$description)])
  dhs_hepg2_avg <- rowMeans(delta[,grepl("HepG2",track_meaning$description) & grepl("DNASE",track_meaning$description)])
  dhs_k562_avg <- rowMeans(delta[,grepl("K562",track_meaning$description) & grepl("DNASE",track_meaning$description)])
  df_summary <- data.frame(cbind(avg,cage_avg,dhs_avg,cage_dhs_avg,cage_hepg2_avg,cage_k562_avg,dhs_hepg2_avg,dhs_k562_avg))
  colnames(df_summary) <- gsub("$",suffix,colnames(df_summary))
  df_summary
}









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




