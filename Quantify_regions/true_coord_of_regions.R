setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep")
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


df <- read.csv("dat_test.csv")
gr <- makeGRangesFromDataFrame(df[,c(2,3,4,5,6)])
gr <- resize(gr,width = 600,fix = "center")
gr$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr)
which(df$seq!=gr$seq)

export.bed(gr,"gr_test.bed")
