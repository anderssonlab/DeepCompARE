setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp")


#CAGEfighR extensions
setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/")
sep <- .Platform$file.sep
cfr.extensions.dir <- "/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/CAGEfightR_extensions"
source(paste(cfr.extensions.dir, "utils.R", sep = sep))
source(paste(cfr.extensions.dir, "support.R", sep = sep))
source(paste(cfr.extensions.dir, "cumulative.R", sep = sep))
source(paste(cfr.extensions.dir, "decompose.R", sep = sep))
source(paste(cfr.extensions.dir, "enhancers.R", sep = sep))
source(paste(cfr.extensions.dir, "noise.R", sep = sep))
source(paste(cfr.extensions.dir, "coverage.R", sep = sep))
source(paste(cfr.extensions.dir, "complexity.R", sep = sep))
source(paste(cfr.extensions.dir, "subsample.R", sep = sep))


#Packages for data manipulation and plotting
library(pheatmap)
library(ggseqlogo)
library(magrittr)
library(ggforce)
library(ggthemes)
library(data.table)
library(assertthat)
library(openxlsx)

#CAGEfighrR and related packages
library(GenomicRanges)
library(SummarizedExperiment)
library(GenomicFeatures)
library(BiocParallel)
library(InteractionSet)
library(Gviz)
library(Rsubread)
library(ggfortify)

#Bioconductor packages for differential expression
library(DESeq2)
library(limma)
library(edgeR)
library(sva)
library(genomation)

#Data packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)






#Rename data packages
bsg <- BSgenome.Hsapiens.UCSC.hg38
txdb<- keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene)

genomeInfo <- keepStandardChromosomes(SeqinfoForUCSCGenome("hg38"))
seqlevels(txdb) <- c(seqlevels(genomeInfo))  #Reorder seqlevels of txdb

register(MulticoreParam(12))





# get CTSSs
cage <- "/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/CAGE/"
mfiles <- list.files(path=cage,pattern="minus.bw",full.names=TRUE)
pfiles <- list.files(path=cage,pattern="plus.bw",full.names=TRUE)
minus.files <- BigWigFileList(mfiles)
plus.files <- BigWigFileList(pfiles)
# Don't use .tpm.bw
names(minus.files) <- as.character(sapply(list.files(path=cage,pattern="minus.bw"), function(n) gsub(".minus.bw","",n)))
names(plus.files) <- as.character(sapply(list.files(path=cage,pattern="plus.bw"), function(n) gsub(".plus.bw","",n)))
CTSSs <- quantifyCTSSs(plusStrand=plus.files,minusStrand=minus.files,genome=genomeInfo)

## Pool re-sequencing runs
n <- as.character(sapply(colnames(CTSSs), function(n) paste(strsplit(n,"_")[[1]][1:2],collapse="_")))
cnt <- assay(CTSSs,"counts")[,match(unique(n),n)]
colnames(cnt) <- unique(n)
for (n in unique(n))
  for (i in which(unique(n)==n)[-1])
    cnt[,n] <- cnt[,n] + assay(CTSSs,"counts")[,i]
merged.CTSSs <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(cnt),
                                                           rowRanges = rowRanges(CTSSs))
SummarizedExperiment::assayNames(merged.CTSSs) <- "counts"
CTSSs <- merged.CTSSs



#Cluster and quantify TCs
CTSSs <- CTSSs %>% calcTPM %>% calcPooled()
TCs <- clusterUnidirectionally(CTSSs)
TCs <- quantifyClusters(CTSSs, TCs)
TCs <- assignTxType(TCs, txModels=txdb)
TCs <- TCs %>% calcTPM()

#Subset for promoters
promoters <- subset(TCs, txType == "promoter")

# calculate cell-type specific expression
expr_mat <- assay(promoters,"counts")
hepg2_avg <- rowMeans(expr_mat[,1:6])
k562_avg <- rowMeans(expr_mat[,7:12])

# remove regions with too low expression
rowData(promoters)$hepg2_expressed <- (hepg2_avg>=10)
rowData(promoters)$k562_expressed <- (k562_avg>=10)
rowData(promoters)$expressed <- rowData(promoters)$hepg2_expressed | rowData(promoters)$k562_expressed
promoters <- subset(promoters,expressed==TRUE)

#Calculate 10-90% IQR
promoters <- calcShape(promoters, CTSSs)

#Look at shape distribution
promoters %>%                                                                   
  rowData %>%                                                                
  as.data.frame %>%                                                          
  ggplot(aes(x=IQR)) +                                                       
  geom_histogram(binwidth=1, fill="orange", alpha=0.75) +                   
  geom_vline(xintercept = 5, linetype="dashed", alpha=0.75, color="black") +
  geom_vline(xintercept = 10, linetype="dashed", alpha=0.75, color="black") +
  facet_zoom(xlim = c(0,100)) +                                              
  labs(x="10-90% IQR", y="Frequency")  

# Define sharp and broad promoters
rowData(promoters)$shape <- ifelse(rowData(promoters)$IQR < 5, "Sharp", "Unknown")
rowData(promoters)$shape[rowData(promoters)$IQR>10] <- "Broad"
table(rowData(promoters)$shape)


# sanity check by plotting seq_logo
promoter_seqs <- promoters %>%                
  rowRanges() %>%                          
  swapRanges() %>%                         
  promoters(upstream=40, downstream=10) %>%
  getSeq(bsg, .)     

promoter_seqs %>%                      
  as.character %>%                   
  split(rowData(promoters)$shape) %>% 
  ggseqlogo(data=., ncol=3, nrow=1) +
  theme_logo() +                     
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())



# get idx of promoters
idx_broad_common <- which(rowData(promoters)$shape=="Broad" & 
                            rowData(promoters)$hepg2_expressed==TRUE & 
                            rowData(promoters)$k562_expressed==TRUE)
idx_sharp_common <- which(rowData(promoters)$shape=="Sharp" & 
                            rowData(promoters)$hepg2_expressed==TRUE & 
                            rowData(promoters)$k562_expressed==TRUE)

idx_broad_hepg2 <- which(rowData(promoters)$shape=="Broad" & 
                          rowData(promoters)$hepg2_expressed==TRUE & 
                          rowData(promoters)$k562_expressed==FALSE)
idx_sharp_hepg2 <- which(rowData(promoters)$shape=="Sharp" & 
                          rowData(promoters)$hepg2_expressed==TRUE & 
                          rowData(promoters)$k562_expressed==FALSE)

idx_broad_k562 <- which(rowData(promoters)$shape=="Broad" & 
                           rowData(promoters)$hepg2_expressed==FALSE & 
                           rowData(promoters)$k562_expressed==TRUE)
idx_sharp_k562 <- which(rowData(promoters)$shape=="Sharp" & 
                           rowData(promoters)$hepg2_expressed==FALSE & 
                           rowData(promoters)$k562_expressed==TRUE)



# subset promoters
promoters_sharp_common <- rowRanges(promoters[idx_sharp_common])
promoters_broad_common <- rowRanges(promoters[idx_broad_common])

promoters_sharp_hepg2 <- rowRanges(promoters[idx_sharp_hepg2])
promoters_broad_hepg2  <- rowRanges(promoters[idx_broad_hepg2])

promoters_sharp_k562 <- rowRanges(promoters[idx_sharp_k562])
promoters_broad_k562 <- rowRanges(promoters[idx_broad_k562])


# resize promoters
promoters_sharp_common <- resize(promoters_sharp_common,fix="center",width=600)
promoters_broad_common <- resize(promoters_broad_common,fix="center",width=600)

promoters_sharp_hepg2 <- resize(promoters_sharp_hepg2,fix="center",width=600)
promoters_broad_hepg2 <- resize(promoters_broad_hepg2,fix="center",width=600)

promoters_sharp_k562 <- resize(promoters_sharp_k562,fix="center",width=600)
promoters_broad_k562 <- resize(promoters_broad_k562,fix="center",width=600)

setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp")




# write to bed file
write.table(as.data.frame(promoters_sharp_common),"promoters_sharp_common.bed",row.names = FALSE,col.names = FALSE)
write.table(as.data.frame(promoters_broad_common),"promoters_broad_common.bed",row.names = FALSE,col.names = FALSE)

write.table(as.data.frame(promoters_sharp_hepg2),"promoters_sharp_hepg2.bed",row.names = FALSE,col.names = FALSE)
write.table(as.data.frame(promoters_broad_hepg2),"promoters_broad_hepg2.bed",row.names = FALSE,col.names = FALSE)

write.table(as.data.frame(promoters_sharp_k562),"promoters_sharp_k562.bed",row.names = FALSE,col.names = FALSE)
write.table(as.data.frame(promoters_broad_k562),"promoters_broad_k562.bed",row.names = FALSE,col.names = FALSE)
