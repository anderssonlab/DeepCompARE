setwd("~/Documents/DeepCompare/SNP_interpretation_raQTL")
source("~/Documents/DeepCompare/Rcode/variant_interpretation.R")

df <- read.csv("Raw_data/hepg2_sure_raQTL.csv")
df$ref.seq <- get_rc_from_char_vec(df$ref.seq)
df$alt.seq <- get_rc_from_char_vec(df$alt.seq)
write.csv(df,"Pd1_reverse_complement_of_raw_data/rc_hepg2_sure_raQTL.csv")

df <- read.csv("Raw_data/k562_sure_raQTL.csv")
df$ref.seq <- get_rc_from_char_vec(df$ref.seq)
df$alt.seq <- get_rc_from_char_vec(df$alt.seq)
write.csv(df,"Pd1_reverse_complement_of_raw_data/rc_k562_sure_raQTL.csv")
