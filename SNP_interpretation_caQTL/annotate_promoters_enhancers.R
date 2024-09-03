library(rtracklayer)
library(CAGEfightR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

df <- read.csv("combined_RASQUAL-results-2PCs_1kbFromPeakCenters-MAF0.1_liver-secondPass-20samples.txt",sep="\t")
df$start <- df$var_position_hg19
df$end <- df$start+1
df <- df[df$var_ID!="SKIPPED",]
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
print(gr)

txdb<- keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genomeInfo <- keepStandardChromosomes(SeqinfoForUCSCGenome("hg19"))
seqlevels(txdb) <- c(seqlevels(genomeInfo))

annot <-  assignTxType(gr, txModels=txdb)
df <- as.data.frame(annot)
write.csv(df,"caQTL_annotated.csv")
