setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions")

library(CAGEfightR)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
source("Rcode/CAGEfightR_extensions/enhancers.R")



seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
genomeInfo <- keepStandardChromosomes(seqInfo,species="Homo_sapiens")

#---------------------------------------
#Step 1: read in all bw files
#---------------------------------------
cage  <- "/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/CAGE"

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

#---------------------------------
# Step 2: split by cell lines
#---------------------------
# param unexpressed: we accept a CTSS only if it has more than UNEXPRESSED counts. unexpressed=1
# param minSamples:we accept a CTSS only if more than MINSAMPLES samples express it. minSamples=0
K562.CTSSs <- subsetBySupport(merged.CTSSs[,grep("K562_",colnames(merged.CTSSs))], inputAssay="counts", unexpressed=1, minSamples=0)
HepG2.CTSSs <- subsetBySupport(merged.CTSSs[,grep("HepG2_",colnames(merged.CTSSs))], inputAssay="counts", unexpressed=1, minSamples=0)


# output cell line CTSSs as RData
CTSSs <- HepG2.CTSSs
colnames(CTSSs)
save(CTSSs,file="/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions/Quant_signal_processed/cage_hepg2.Rdata")
CTSSs <- K562.CTSSs
colnames(CTSSs)
save(CTSSs,file="/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions/Quant_signal_processed/cage_k562.Rdata")


#----------------------------------------------------------
# Step 3: calculate divergentLoci,
# output resulting region for reproducibility analysis
#--------------------------------------------------------

export_RE <- function(CTSSs,expr_thresh,support_thresh,out_name){
  CTSSs <- calcTPM(CTSSs, inputAssay="counts", outputAssay="TPM")
  CTSSs <- calcPooled(CTSSs,inputAssay="TPM")
  TCs<- clusterUnidirectionally(CTSSs)
  re <- divergentLoci(TCs,CTSSs) 
  
  expr <- quantifyDivergentLoci(re, CTSSs, requireDisjoint=FALSE)
  expressed_above_noise <- assay(expr,"counts")>=expr_thresh
  idx <- which(rowSums(expressed_above_noise)>=support_thresh)
  
  print(length(idx))
  gr <- rowRanges(expr)[idx]
  mcols(gr) <- NULL
  print(length(gr))
  export.bed(gr,out_name)
}

export_RE(HepG2.CTSSs,6,2,"CAGE/Peaks_reproducible_within_modality/CAGE_HepG2.bed")
export_RE(K562.CTSSs,6,2,"CAGE/Peaks_reproducible_within_modality/CAGE_K562.bed")

export_RE(HepG2.CTSSs,6,1,"CAGE/Peaks_irreproducible/CAGE_HepG2_irreproducible.bed")
export_RE(K562.CTSSs,6,1,"CAGE/Peaks_irreproducible/CAGE_K562_irreproducible.bed")

# Here loose contains reproducible and irreproducible by definition
export_RE(HepG2.CTSSs,1,1,"/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/CAGE_HepG2_loose.bed")
export_RE(K562.CTSSs,1,1,"/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/CAGE_K562_loose.bed")







#----------------------------------------------------------
# Archived
#----------------------------------------------------------


# Step 3, alternative: quantify on cell line DHS 


genomeInfo <- keepStandardChromosomes(seqInfo,species="Homo_sapiens")
negative <- import("/isdata/alab/people/pcr980/Raw_data/Bed_processed/hg38_non_hotspots_non_blacklist_600bp_windows_intergenic.bed", format="BED", genome=genomeInfo)
dhs_hepg2 <- import.bed("/isdata/alab/people/pcr980/Raw_data/Bed_processed/DHS_HepG2.bed",genome=genomeInfo)
dhs_k562 <- import.bed("/isdata/alab/people/pcr980/Raw_data/Bed_processed/DHS_K562.bed",genome=genomeInfo)

negative$thick <- IRanges(start=ceiling((start(negative)+end(negative))/2), width=1)
dhs_hepg2$thick <- IRanges(start=ceiling((start(dhs_hepg2)+end(dhs_hepg2))/2), width=1)
dhs_k562$thick <- IRanges(start=ceiling((start(dhs_k562)+end(dhs_k562))/2), width=1)

dhs_hepg2 <- dhs_hepg2[1]
strand(dhs_hepg2) <- "*"
expr_dl <- quantifyDivergentLoci(dhs_hepg2, HepG2.CTSSs, requireDisjoint=FALSE)
assay(expr_dl,"counts")

expr_cluster <- quantifyClusters(HepG2.CTSSs,clusters = dhs_hepg2)
assay(expr_cluster,"counts")



eHepG2.DHS.expr <- quantifyDivergentLoci(dhs_hepg2, HepG2.CTSSs, requireDisjoint=FALSE)
HepG2.negative.expr <- quantifyDivergentLoci(negative, HepG2.CTSSs, requireDisjoint=FALSE)
HepG2.thres <- apply(assay(HepG2.negative.expr,"counts"),2,function(x) quantile(x,probs=0.99))
names(HepG2.thres) <- colnames(HepG2.negative.expr)


K562.DHS.expr <- quantifyDivergentLoci(dhs_k562, K562.CTSSs, requireDisjoint=FALSE)
K562.negative.expr <- quantifyDivergentLoci(negative, K562.CTSSs, requireDisjoint=FALSE)
K562.thres <- apply(assay(K562.negative.expr,"counts"),2,function(x) quantile(x,probs=0.99))
names(K562.thres) <- colnames(K562.negative.expr)


num.rep <- 2
assay(HepG2.DHS.expr,"expressed") <- Matrix::Matrix(sapply(names(HepG2.thres), function(n) assay(HepG2.DHS.expr,"counts")[,n] > HepG2.thres[n]),sparse=FALSE)
rowData(HepG2.DHS.expr)[,"N_expressed"] <- rowSums(assay(HepG2.DHS.expr,"expressed")[,grep("_N",names(HepG2.thres))])>=num.rep
assay(K562.DHS.expr,"expressed") <- Matrix::Matrix(sapply(names(K562.thres), function(n) assay(K562.DHS.expr,"counts")[,n] > K562.thres[n]),sparse=FALSE)
rowData(K562.DHS.expr)[,"N_expressed"] <- rowSums(assay(K562.DHS.expr,"expressed")[,grep("_N",names(K562.thres))])>=num.rep

export.bed(rowRanges(HepG2.DHS.expr)[rowRanges(HepG2.DHS.expr)$N_expressed],
           "/isdata/alab/people/pcr980/Raw_data_CAGE/QDHS_CAGE_HepG2.bed")
export.bed(rowRanges(K562.DHS.expr)[rowRanges(K562.DHS.expr)$N_expressed],
           "/isdata/alab/people/pcr980/Raw_data_CAGE/QDHS_CAGE_K562.bed")




# calculate TPM expression
K562.DHS.expr <- calcTPM(K562.DHS.expr)
rowData(K562.DHS.expr)[,"N_expression"] <- rowMeans(assay(K562.DHS.expr,"TPM")[,grep("_N",names(K562.thres))])

HepG2.DHS.expr <- calcTPM(HepG2.DHS.expr)
rowData(HepG2.DHS.expr)[,"N_expression"] <- rowMeans(assay(HepG2.DHS.expr,"TPM")[,grep("_N",names(HepG2.thres))])




