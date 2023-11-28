setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/Rcode/Function_utils.R")
library(readr)
library(liftOver)

data_dir <- "SuRE/NarrowPeak_hg19_c4/"
EDA_dir <- "EDAs/"
black_dir <- "SuRE_files/SuRE_SNP/"
ch <- import.chain("/maps/projects/ralab/people/pcr980/Resource_Genome_annot/hg19ToHg38.over.chain")
#-----------------------------------------------------------
# EDA on each file
# Result: found no difference in distribution of peak length
#-----------------------------------------------------------

pdf(paste0(EDA_dir,"0.SuRE_MACS3_output.pdf"),height = 7, width = 22)
lapply(list.files(data_dir,full.name=T), function(file_path){
  file_name <- full_path_to_file_name(file_path)
  message(paste0("EDA on ",file_name))
  gr <- read_file(file_path,remove_mcols=F)
  par(mfrow=c(1,3))
  hist(log2(width(gr)),
       main=paste0(file_name,":log2 width"),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(width(gr)[width(gr)<2000],
       main=paste0(file_name,":width zoom in"),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(log2(calculate_distance_self(gr)),
       main=paste0(file_name,":Dist2Nearest"),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
})
dev.off()

#-------------------------------------------------------------------
# remove SNP, remove excessively long peaks
# Result dir: NarrowPeak_hg19_RMSNP/ (to deduce regions reproducibile within modality)
#------------------------------------------------------------------


out_dir <- "SuRE/NarrowPeak_hg19_c4_RMSNP/"

lapply(list.files(data_dir), function(file){
  gr <- read_file(paste0(data_dir,file),remove_mcols=F)
  gr <- gr[width(gr)<2000]
  # remove regions overlapping with SNP
  black_file <- list.files(black_dir,pattern=paste0(parse_number(file),"_hg19.bed"))
  black <- unique(import(paste0(black_dir,black_file),format="bed"))
  ov <- findOverlaps(gr,black,ignore.strand=TRUE)
  gr <- gr[-unique(queryHits(ov))]
  write_narrowPeak(gr,paste0(out_dir,file))
})

#----------------------------------------------------------------------------------------------
# Lift hg19 to hg38, remove SNP 
# Result dir: NarrowPeak_hg38_RMSNP/ (to rescue irreproducible sequences of other modalities)
#----------------------------------------------------------------------------------------------

out_dir <- "SuRE/NarrowPeak_hg38_c4_RMSNP/"

lapply(list.files(data_dir), function(file){
  gr <- read_file(paste0(data_dir,file),remove_mcols=T)
  # remove regions overlapping with SNP
  black_file <- list.files(black_dir,pattern=paste0(parse_number(file),"_hg19.bed"))
  black <- unique(import(paste0(black_dir,black_file),format="bed"))
  ov <- findOverlaps(gr,black,ignore.strand=TRUE)
  gr <- gr[-unique(queryHits(ov))]
  # lift over to hg38
  gr <- unlist(liftOver(gr, ch))
  out_name <- gsub("hg19","hg38",file)
  write_narrowPeak(gr,paste0(out_dir,out_name))
})



#----------------------------------------------------------------------------------------------
# Lift hg19 to hg38, don't remove SNP 
# Result dir: NarrowPeak_hg38_c2/ (to deduce loose)
#----------------------------------------------------------------------------------------------
data_dir <- "SuRE/NarrowPeak_hg19_c2/"
out_dir <- "SuRE/NarrowPeak_hg38_c2/"

lapply(list.files(data_dir), function(file){
  gr <- read_file(paste0(data_dir,file),remove_mcols=T)
  gr <- unlist(liftOver(gr, ch))
  out_name <- gsub("hg19","hg38",file)
  write_narrowPeak(gr,paste0(out_dir,out_name))
})

