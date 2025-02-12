setwd("~/Documents/DeepCompare/Fig5_DNase")


## Load library and data
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)
library(readr)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters <- promoters(txdb,upstream=500,downstream=500)
cds <- unlist(cdsBy(txdb,"gene"))



# specify ROADMAP file name and output suffix
FILE_NAME <- "DNase_E123_k562.bed"
OUT_SUFFIX <- "k562"



# read data
gr_dnase <- import(FILE_NAME,format="bed")
# resize to 600bp
gr_dnase <- resize(gr_dnase,width=600,fix="center")
# remove regions overlapping with cds
ov <- findOverlaps(gr_dnase,cds,minoverlap=1,ignore.strand=TRUE)
gr_dnase <- gr_dnase[-queryHits(ov)]


# annotate proximal / distal
ov <- findOverlaps(gr_dnase,promoters,minoverlap=1,ignore.strand=TRUE)
gr_dnase$location <- "distal"
gr_dnase[queryHits(ov)]$location <- "proximal"



# annotate tissue invariant
dhs_idx <- read_tsv("DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
dhs_ti <- dhs_idx[dhs_idx$component == "Tissue invariant", ]
dhs_ti <- GRanges(seqnames = dhs_ti$seqname, ranges = IRanges(start = dhs_ti$start, end = dhs_ti$end))
ov <- findOverlaps(gr_dnase,dhs_ti,ignore.strand=TRUE)
gr_dnase$tissue_invariance <- "no"
gr_dnase[queryHits(ov)]$tissue_invariance <- "yes"


# calculate constraints
constraints <- read.table(gzfile("./constraint_z_genome_1kb.qc.download.txt.gz"), header = TRUE)
constraints <- GRanges(constraints)

# split dhs by constraints
calc_constraints<- function(dhs) {
  # Find overlaps
  ov <- findOverlaps(dhs, constraints)
  
  # Calculate overlap lengths
  overlap_lengths <-width(pintersect(dhs[queryHits(ov)], constraints[subjectHits(ov)]))
  df_ov <- data.frame(queryHits(ov), subjectHits(ov), overlap_lengths)
  names(df_ov) <- c("dhs_idx", "constraint_idx", "overlap_length")
  
  # Extract the corresponding z scores and weights
  df_ov$z_scores <- constraints$z[df_ov$constraint_idx]
  df_ov$weighted_z <- df_ov$z_scores * df_ov$overlap_length / 600
  
  # group by dhs_idx, sum the weighted z scores
  weighted_z <- tapply(df_ov$weighted_z, df_ov$dhs_idx, sum)
  
  # Add the weighted z-scores back to dhs
  mcols(dhs)$weighted_z <- NA
  mcols(dhs)$weighted_z[as.numeric(names(weighted_z))] <- weighted_z
  
  return(dhs)
}



gr_dnase <- calc_constraints(gr_dnase)
# retain only chr1-22,X,Y
df_dnase <- as.data.frame(gr_dnase)
df_dnase <- df_dnase[df_dnase$seqnames %in% paste0("chr",c(1:22,"X","Y")),]

df_proximal <- df_dnase[df_dnase$location == "proximal",]
df_distal <- df_dnase[df_dnase$location == "distal",]


write.table(df_proximal,paste0("dhs_proximal_",OUT_SUFFIX,".tsv"),row.names=FALSE)
write.table(df_distal,paste0("dhs_distal_",OUT_SUFFIX,".tsv"),row.names=FALSE)




