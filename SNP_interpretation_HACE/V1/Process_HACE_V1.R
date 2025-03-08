setwd("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_HACE/V1")
library(readxl)
library(Biostrings)
df <- as.data.frame(read_excel("Log2ORCD69hi-to-low.xlsx"))
ref_seq <- DNAString(paste0(df[,1],collapse = ""))
reletive_position <- df[,2]-9764855           
c_idx <- which(!is.na(df[,3]))      
g_idx <- which(!is.na(df[,5]))
seqs_c2t <- lapply(c_idx,function(this_idx){
  replaceLetterAt(ref_seq,this_idx,"T")
})

seqs_c2t <- DNAStringSet(seqs_c2t)
logratio_c2t_1 <- df[c_idx,3]
logratio_c2t_2 <- df[c_idx,4]
relevant_position_c2t <- reletive_position[c_idx]

seqs_g2a <- lapply(g_idx,function(this_idx){
  replaceLetterAt(ref_seq,this_idx,"A")
})
seqs_g2a <- DNAStringSet(seqs_g2a)
logratio_g2a_1 <- df[g_idx,5]
logratio_g2a_2 <- df[g_idx,6]
relevant_position_g2a <- reletive_position[g_idx]



seqs_df <- as.data.frame(c(seqs_c2t,seqs_g2a))
seqs_df$log2ratio_1 <- c(logratio_c2t_1,logratio_g2a_1)
seqs_df$log2ratio_2 <- c(logratio_c2t_2,logratio_g2a_2)
seqs_df$relevant_position <- c(relevant_position_c2t,relevant_position_g2a)
colnames(seqs_df)[1] <- "alt_seq"
seqs_df$ref_seq <- as.character(ref_seq)


write.csv(seqs_df,"Validation_SNP_HACE/log2ratio_clean.csv",row.names = F)

df_pred <- read.csv("Validation_SNP_HACE/prediction_log2ratio.csv")

plot(df_pred$log1p_signal_sure_k562_alt-df_pred$log1p_signal_sure_k562_ref,seqs_df$log2ratio_2)
