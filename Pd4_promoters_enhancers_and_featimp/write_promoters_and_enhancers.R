library(rtracklayer)
library(CAGEfightR)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd("/maps/projects/ralab/people/pcr980/DeepCompare")


#----------------------------------------
# Promoters: Annotate CAGE peaks with txdb
#----------------------------------------
txdb<- keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genomeInfo <- keepStandardChromosomes(SeqinfoForUCSCGenome("hg38"))
seqlevels(txdb) <- c(seqlevels(genomeInfo))


bed <- import.bed("Pd1_bed_processed/CAGE_HepG2.bed")
annot <-  assignTxType(bed, txModels=txdb)
promoters <- annot[annot$txType=="promoter"]
promoters <- resize(promoters,width=600,fix="center")
export.bed(unique(promoters),"Pd4_promoters_enhancers_and_featimp/promoters_hepg2.bed")


bed <- import.bed("Pd1_bed_processed/CAGE_K562.bed")
annot <-  assignTxType(bed, txModels=txdb)
promoters <- annot[annot$txType=="promoter"]
promoters <- resize(promoters,width=600,fix="center")
export.bed(unique(promoters),"Pd4_promoters_enhancers_and_featimp/promoters_k562.bed")



#--------------------------------------------------------------------
# Enhancers: subset by DHS, check histone marks (H3K4me1) before and after
#--------------------------------------------------------------------

get_enhancer_from_dhs_and_starr(dhs_file){
  dhs <- import.bed("Pd1_bed_processed/DHS_K562.bed")
  names(dhs) <- paste0("Seq",1:length(dhs))
  starr <- import.bed("Pd1_bed_processed/STARR_K562.bed")
  names(starr) <- paste0("Seq",1:length(starr))
  starr <- subsetByOverlaps(starr,dhs)
  starr <- resize(starr,fix="center",width=600)
  export.bed(unique(starr),"Pd4_promoters_enhancers_and_featimp/enhancers_k562.bed")
}

get_enhancer_from_dhs_and_starr("Pd1_bed_processed/DHS_K562.bed",)






#---------------------
# Archived: check histone signal of selected regions.
#---------------------

source("Scripts_R/Generate_data/Function_bigwig.R")
bw_path <- "Raw_data/Histone/H3K4me1_K562.bigWig"
dhs_histone <- quantify_one_bw_file(dhs,bw_path)
starr_histone <- quantify_one_bw_file(starr,bw_path)
data <- data.frame(
  value = c(dhs_histone[,1], starr_histone[,1], enhancer_histone[,1]),
  group = c(rep("DHS", length(dhs_histone[,1])),
            rep("STARR", length(starr_histone[,1])),
            rep("Enhancer", length(enhancer_histone[,1])))
)

# Plot the density plot
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Value", y = "Density", title = "Density Plot of DHS, STARR, and Enhancer Histones") +
  scale_fill_manual(values = c("DHS" = "blue", "STARR" = "red", "Enhancer" = "green"))

