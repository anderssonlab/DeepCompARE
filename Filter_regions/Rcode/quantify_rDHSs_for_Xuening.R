
library(CAGEfightR)
source("CAGEfightR_extensions/utils.R")
library(rtracklayer)
library(GenomicRanges)
#library(BSgenome.Hsapiens.UCSC.hg38)


    ##seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
    ##save(seqInfo,file="seqInfo.RData")
    load("seqInfo.RData")
genomeInfo <- keepStandardChromosomes(seqInfo,species="Homo_sapiens")
negative <- import("/isdata/alab/people/masa/PROCESSED_PEAKS/DNAse-Encode-rDHS/original/hg38_non_hotspots_non_blacklist_600bp_windows_intergenic.bed", format="BED", genome=genomeInfo)
negative$thick <- IRanges(start=ceiling((start(negative)+end(negative))/2), width=1)

rDHSs <- import("/isdata/alab/people/pcr980/Raw_bed_files/rDHS_600bp.bed", format="BED", genome=genomeInfo)
##rDHSs <- import("/isdata/alab/people/masa/PROCESSED_PEAKS/DNAse-Encode-rDHS/original/V2.hg38-rDHS-Filtered-region.bed", format="BED", genome=genomeInfo)
rDHSs$thick <- IRanges(start=ceiling((start(rDHSs)+end(rDHSs))/2), width=1)

## Refocus to do CAGE quantification in the 300bp central regions
start(rDHSs) <- start(rDHSs$thick)-150
end(rDHSs) <- start(rDHSs$thick)+150

start(negative) <- start(negative$thick)-150
end(negative) <- start(negative$thick)+150

## Read in CAGE data
cage  <- "/isdata/alab/data/CAGE/cellLines/bw_files/"

mfiles <- list.files(path=cage,pattern="minus.bw",full.names=TRUE)
pfiles <- list.files(path=cage,pattern="plus.bw",full.names=TRUE)
minus.files <- BigWigFileList(mfiles)
plus.files <- BigWigFileList(pfiles)
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

## Split in cell line CTSSs
K562.CTSSs <- subsetBySupport(merged.CTSSs[,grep("K562",colnames(merged.CTSSs))], inputAssay="counts", unexpressed=0, minSamples=0)
HCT116.CTSSs <- subsetBySupport(merged.CTSSs[,grep("HCT116",colnames(merged.CTSSs))], inputAssay="counts", unexpressed=0, minSamples=0)
A549.CTSSs <- subsetBySupport(merged.CTSSs[,grep("A549",colnames(merged.CTSSs))], inputAssay="counts", unexpressed=0, minSamples=0)
HepG2.CTSSs <- subsetBySupport(merged.CTSSs[,grep("HepG2",colnames(merged.CTSSs))], inputAssay="counts", unexpressed=0, minSamples=0)
GM12878.CTSSs <- readRDS("GM12878/20210526_ctss_N_1mio.rds")
seqlevels(GM12878.CTSSs,pruning.mode="coarse") <- seqlevels(genomeInfo)

K562.rDHSs.expr <- quantifyDivergentLoci(rDHSs, K562.CTSSs, requireDisjoint=FALSE)
K562.negative.expr <- quantifyDivergentLoci(negative, K562.CTSSs, requireDisjoint=FALSE)

HCT116.rDHSs.expr <- quantifyDivergentLoci(rDHSs, HCT116.CTSSs, requireDisjoint=FALSE)
HCT116.negative.expr <- quantifyDivergentLoci(negative, HCT116.CTSSs, requireDisjoint=FALSE)

A549.rDHSs.expr <- quantifyDivergentLoci(rDHSs, A549.CTSSs, requireDisjoint=FALSE)
A549.negative.expr <- quantifyDivergentLoci(negative, A549.CTSSs, requireDisjoint=FALSE)

HepG2.rDHSs.expr <- quantifyDivergentLoci(rDHSs, HepG2.CTSSs, requireDisjoint=FALSE)
HepG2.negative.expr <- quantifyDivergentLoci(negative, HepG2.CTSSs, requireDisjoint=FALSE)

GM12878.rDHSs.expr <- quantifyDivergentLoci(rDHSs, GM12878.CTSSs, requireDisjoint=FALSE)
GM12878.negative.expr <- quantifyDivergentLoci(negative, GM12878.CTSSs, requireDisjoint=FALSE)


## Identify noise threshold
K562.thres <- apply(assay(K562.negative.expr,"counts"),2,function(x) quantile(x,probs=0.999))
names(K562.thres) <- colnames(K562.negative.expr)

HepG2.thres <- apply(assay(HepG2.negative.expr,"counts"),2,function(x) quantile(x,probs=0.999))
names(HepG2.thres) <- colnames(HepG2.negative.expr)

A549.thres <- apply(assay(A549.negative.expr,"counts"),2,function(x) quantile(x,probs=0.999))
names(A549.thres) <- colnames(A549.negative.expr)

HCT116.thres <- apply(assay(HCT116.negative.expr,"counts"),2,function(x) quantile(x,probs=0.999))
names(HCT116.thres) <- colnames(HCT116.negative.expr)

GM12878.thres <- apply(assay(GM12878.negative.expr,"counts"),2,function(x) quantile(x,probs=0.999))
names(GM12878.thres) <- colnames(GM12878.negative.expr)


## Binarize replicate expression by noise threshold
num.rep <- 1
assay(K562.rDHSs.expr,"expressed") <- Matrix::Matrix(sapply(names(K562.thres), function(n) assay(K562.rDHSs.expr,"counts")[,n] > K562.thres[n]),sparse=FALSE)
rowData(K562.rDHSs.expr)[,"N_expressed"] <- rowSums(assay(K562.rDHSs.expr,"expressed")[,grep("_N",names(K562.thres))])>=num.rep

assay(HepG2.rDHSs.expr,"expressed") <- Matrix::Matrix(sapply(names(HepG2.thres), function(n) assay(HepG2.rDHSs.expr,"counts")[,n] > HepG2.thres[n]),sparse=FALSE)
rowData(HepG2.rDHSs.expr)[,"N_expressed"] <- rowSums(assay(HepG2.rDHSs.expr,"expressed")[,grep("_N",names(HepG2.thres))])>=num.rep

assay(HCT116.rDHSs.expr,"expressed") <- Matrix::Matrix(sapply(names(HCT116.thres), function(n) assay(HCT116.rDHSs.expr,"counts")[,n] > HCT116.thres[n]),sparse=FALSE)
rowData(HCT116.rDHSs.expr)[,"N_expressed"] <- rowSums(assay(HCT116.rDHSs.expr,"expressed")[,grep("_N",names(HCT116.thres))])>=num.rep

assay(A549.rDHSs.expr,"expressed") <- Matrix::Matrix(sapply(names(A549.thres), function(n) assay(A549.rDHSs.expr,"counts")[,n] > A549.thres[n]),sparse=FALSE)
rowData(A549.rDHSs.expr)[,"N_expressed"] <- rowSums(assay(A549.rDHSs.expr,"expressed")[,grep("_N",names(A549.thres))])>=num.rep

assay(GM12878.rDHSs.expr,"expressed") <- Matrix::Matrix(sapply(names(GM12878.thres), function(n) assay(GM12878.rDHSs.expr,"counts")[,n] > GM12878.thres[n]),sparse=FALSE)
rowData(GM12878.rDHSs.expr)[,"N_expressed"] <- rowSums(assay(GM12878.rDHSs.expr,"expressed")[,grep("N_",names(GM12878.thres))])>=num.rep


## Pool replicates and average TPM values
K562.rDHSs.expr <- calcTPM(K562.rDHSs.expr)
rowData(K562.rDHSs.expr)[,"N_expression"] <- rowMeans(assay(K562.rDHSs.expr,"TPM")[,grep("_N",names(K562.thres))])

HepG2.rDHSs.expr <- calcTPM(HepG2.rDHSs.expr)
rowData(HepG2.rDHSs.expr)[,"N_expression"] <- rowMeans(assay(HepG2.rDHSs.expr,"TPM")[,grep("_N",names(HepG2.thres))])

A549.rDHSs.expr <- calcTPM(A549.rDHSs.expr)
rowData(A549.rDHSs.expr)[,"N_expression"] <- rowMeans(assay(A549.rDHSs.expr,"TPM")[,grep("_N",names(A549.thres))])

HCT116.rDHSs.expr <- calcTPM(HCT116.rDHSs.expr)
rowData(HCT116.rDHSs.expr)[,"N_expression"] <- rowMeans(assay(HCT116.rDHSs.expr,"TPM")[,grep("_N",names(HCT116.thres))])

GM12878.rDHSs.expr <- calcTPM(GM12878.rDHSs.expr)
rowData(GM12878.rDHSs.expr)[,"N_expression"] <- rowMeans(assay(GM12878.rDHSs.expr,"TPM")[,grep("N_",names(GM12878.thres))])


## Overall positive and negative acording to expression

rowData(K562.rDHSs.expr)[,"negative"] <- rowSums(assay(K562.rDHSs.expr,"expressed"))==0
rowData(A549.rDHSs.expr)[,"negative"] <- rowSums(assay(A549.rDHSs.expr,"expressed"))==0
rowData(HCT116.rDHSs.expr)[,"negative"] <- rowSums(assay(HCT116.rDHSs.expr,"expressed"))==0
rowData(HepG2.rDHSs.expr)[,"negative"] <- rowSums(assay(HepG2.rDHSs.expr,"expressed"))==0
rowData(GM12878.rDHSs.expr)[,"negative"] <- rowSums(assay(GM12878.rDHSs.expr,"expressed"))==0

rowData(K562.rDHSs.expr)[,"positive"] <- rowData(K562.rDHSs.expr)[,"N_expressed"]
rowData(A549.rDHSs.expr)[,"positive"] <- rowData(A549.rDHSs.expr)[,"N_expressed"]
rowData(HCT116.rDHSs.expr)[,"positive"] <- rowData(HCT116.rDHSs.expr)[,"N_expressed"]
rowData(HepG2.rDHSs.expr)[,"positive"] <- rowData(HepG2.rDHSs.expr)[,"N_expressed"]
rowData(GM12878.rDHSs.expr)[,"positive"] <- rowData(GM12878.rDHSs.expr)[,"N_expressed"]

rDHS.class <- do.call("cbind",
                      list(
                          rowData(K562.rDHSs.expr)[,"positive"],
                          rowData(A549.rDHSs.expr)[,"positive"],
                          rowData(HCT116.rDHSs.expr)[,"positive"],
                          rowData(HepG2.rDHSs.expr)[,"positive"],
                          rowData(GM12878.rDHSs.expr)[,"positive"]))

colnames(rDHS.class) <- c("K562","A549","HCT116","HepG2","GM12878")
rownames(rDHS.class) <- rowData(K562.rDHSs.expr)$name

save(rDHS.class,file="rDHS_data/rDHS.CAGE.class.RData")
write.csv(as.data.frame(rDHS.class),file="rDHS_data/rDHS.CAGE.class.csv")

## Save data and write expression to file
save(K562.CTSSs,K562.negative.expr,K562.thres,K562.rDHSs.expr,file="rDHS_data/K562_data_CAGE.RData")
write.csv(as.data.frame(rowRanges(K562.rDHSs.expr)),file="rDHS_data/K562_rDHSs_data_CAGE.csv")

save(HCT116.CTSSs,HCT116.negative.expr,HCT116.thres,HCT116.rDHSs.expr,file="rDHS_data/HCT116_data_CAGE.RData")
write.csv(as.data.frame(rowRanges(HCT116.rDHSs.expr)),file="rDHS_data/HCT116_rDHSs_data_CAGE.csv")

save(HepG2.CTSSs,HepG2.negative.expr,HepG2.thres,HepG2.rDHSs.expr,file="rDHS_data/HepG2_data_CAGE.RData")
write.csv(as.data.frame(rowRanges(HepG2.rDHSs.expr)),file="rDHS_data/HepG2_rDHSs_data_CAGE.csv")

save(A549.CTSSs,A549.negative.expr,A549.thres,A549.rDHSs.expr,file="rDHS_data/A549_data_CAGE.RData")
write.csv(as.data.frame(rowRanges(A549.rDHSs.expr)),file="rDHS_data/A549_rDHSs_data_CAGE.csv")

save(GM12878.CTSSs,GM12878.negative.expr,GM12878.thres,GM12878.rDHSs.expr,file="rDHS_data/GM12878_data_CAGE.RData")
write.csv(as.data.frame(rowRanges(GM12878.rDHSs.expr)),file="rDHS_data/GM12878_rDHSs_data_CAGE.csv")

