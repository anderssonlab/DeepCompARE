setwd("/maps/projects/ralab/people/pcr980/DeepCompare/MPRA_combinatorics")

df <- read.csv("Raw_data/my_test_MPRA.csv")
sequences <- c(df$ref_seq,df$alt_seq)
idx <- paste0("Seq",1:length(sequences))
df <- as.data.frame(cbind(idx,sequences))

write.csv(df,"Pd1_formatted_sequences/MPRA_seqs.csv")
