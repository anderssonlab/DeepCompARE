library(rtracklayer)
library(GenomicRanges)

# Read the BED file into a GRanges object
bed_disease <- import("ClinVar_non_benign_non_coding.bed", format = "BED")
bed_disease <- import("whole_blood_eQTLs_PIP_05_non_coding.bed", format = "BED")
bed_disease <- import("red_blood_cell_trait_GWAS_PIP_02_non_coding.bed", format = "BED")
bed_disease <- import("background_SNPs_non_coding.bed", format = "BED")

# read annotations
enhancers <- import("enhancers_k562.bed", format = "BED")
promoters <- import("promoters_k562.bed", format = "BED")

codependent <- read.csv("tfbs_k562_codependent.csv")
codependent <- GRanges(seqnames = codependent$chromosome, ranges = IRanges(start = codependent$start, end = codependent$end))

redundant <- read.csv("tfbs_k562_redundant.csv")
redundant <- GRanges(seqnames = redundant$chromosome, ranges = IRanges(start = redundant$start, end = redundant$end))


# calculate overlaps
ov <- findOverlaps(bed_disease, enhancers)
bed_disease$on_enhancer <- FALSE
bed_disease$on_enhancer[queryHits(ov)] <- TRUE

ov <- findOverlaps(bed_disease, promoters)
bed_disease$on_promoter <- FALSE
bed_disease$on_promoter[queryHits(ov)] <- TRUE

ov <- findOverlaps(bed_disease, codependent)
bed_disease$on_codependent <- FALSE
bed_disease$on_codependent[queryHits(ov)] <- TRUE

ov <- findOverlaps(bed_disease, redundant)
bed_disease$on_redundant <- FALSE
bed_disease$on_redundant[queryHits(ov)] <- TRUE


# count 
enhancer_redundant=sum(bed_disease$on_enhancer & bed_disease$on_redundant)
enhancer_codependent=sum(bed_disease$on_enhancer & bed_disease$on_codependent)

promoter_redundant=sum(bed_disease$on_promoter & bed_disease$on_redundant)
promoter_codependent=sum(bed_disease$on_promoter & bed_disease$on_codependent)




