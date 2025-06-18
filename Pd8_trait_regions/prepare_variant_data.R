library(GenomicRanges)
library(readr)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

SNP.file <- "../nucCAGE/0.External_resources/all.bg.SNPs.hg38.baseline.v1.1.bed.sorted"
eQTL.file <- "../nucCAGE/0.External_resources/GTEx_30tissues_release1_pip_01_Whole_Blood.tsv"
GWAS.file <- "../nucCAGE/0.External_resources/UKBB_94traits_release1_pip01_clean_uniq_hg38.bed"
ClinVar.file <- "../nucCAGE/0.External_resources/ClinVar_GRCh38_nonbenign_noncoding_SNV.tsv"

eQTL.pip.thres <- 0.5
GWAS.pip.thres <- 0.2
traits <- c("HbA1c","Hb","RBC","MCV","MCH","MCHC")

## Annotated CDS
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cds <- unlist(cdsBy(txdb,"gene"))

## Background SNPs (9,913,974 unique autosomal SNP locations)
SNP.data <- import.bed(SNP.file)
SNP.data <- sort(SNP.data)
olap <- findOverlaps(SNP.data,cds,minoverlap=1)
SNP.data$name <- NULL
SNP.data <- unique(SNP.data[-queryHits(olap)])

## eQTL data
eQTL.data <- read_tsv(eQTL.file,col_names=FALSE,show_col_types=FALSE)[,c(1,3)]
colnames(eQTL.data) <- c("variant","pip")
pos <- as.numeric(sapply(eQTL.data$variant, function(n) strsplit(n,"_")[[1]][2]))+1
eQTL.data <- GRanges(seqnames=as.character(sapply(eQTL.data$variant, function(n) strsplit(n,"_")[[1]][1])),
                     score=eQTL.data$pip,
                     ranges=IRanges(pos, pos))
eQTL.data <- sort(eQTL.data)
seqlevels(eQTL.data,pruning.mode="coarse") <- seqlevels(SNP.data)
eQTL.data <- subset(eQTL.data, score >= eQTL.pip.thres)
eQTL.data$score <- NULL
olap <- findOverlaps(eQTL.data,cds,minoverlap=1)
eQTL.data <- unique(eQTL.data[-queryHits(olap)])

## GWAS data
GWAS.data <- read_tsv(GWAS.file,col_names=FALSE,show_col_types=FALSE)
colnames(GWAS.data) <- c("chr","start","end","name","score","strand","RS_id","trait","pip")
GWAS.data <- makeGRangesFromDataFrame(GWAS.data,keep.extra.columns = TRUE, starts.in.df.are.0based=TRUE)
GWAS.data <- subset(GWAS.data, pip >= GWAS.pip.thres)
seqlevels(GWAS.data,pruning.mode="coarse") <- seqlevels(SNP.data)
olap <- findOverlaps(GWAS.data,cds,minoverlap=1)
GWAS.data <- GWAS.data[-queryHits(olap)]
GWAS.data <- subset(GWAS.data, trait %in% traits)
GWAS.data$pip <- NULL
GWAS.data$trait <- NULL
GWAS.data <- unique(GWAS.data)

## ClinVar data
ClinVar.data <- read_tsv(ClinVar.file,show_col_types=FALSE)
ClinVar.data <- makeGRangesFromDataFrame(ClinVar.data)
seqlevels(ClinVar.data,pruning.mode="coarse") <- seqlevels(SNP.data)
olap <- findOverlaps(ClinVar.data,cds,minoverlap=1)
ClinVar.data <- ClinVar.data[-queryHits(olap)]
ClinVar.data <- sort(ClinVar.data)
ClinVar.data <- unique(ClinVar.data)

## Write to bed files
export.bed(SNP.data,"background_SNPs_non_coding.bed")
export.bed(eQTL.data,"whole_blood_eQTLs_PIP_05_non_coding.bed")
export.bed(GWAS.data,"red_blood_cell_trait_GWAS_PIP_02_non_coding.bed")
export.bed(ClinVar.data,"ClinVar_non_benign_non_coding.bed")
