# Basic functions for reading in files
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(dplyr)
library(stringr)

read_file <- function(path,remove_mcols=FALSE){
  # read in BED, narrowPeak, and output of STARRPeaker, output GRanges
  if (grepl("bed",path)){
    file <- unique(import(path,format="bed"))
  }
  if (grepl("narrowPeak",path)){
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    file <- unique(import(path,format="BED",extraCols = extraCols_narrowPeak))
  }
  if (grepl("starrPeak",path)){
    file <- read.table(path)
    if (dim(file)[2]==11){
      colnames(file) <- c("chr","start","end","name","score","strand","fc","input_coverage","output_coverage","pval","qval")
      file <- file[,c("chr","start","end","strand","name","score","fc","input_coverage","output_coverage","pval","qval")]
      file <- unique(makeGRangesFromDataFrame(file,keep.extra.columns=T))
    }
    if (dim(file)[2]==10){
      colnames(file) <- c("chr","start","end","name","score","strand","fc","output_coverage","pval","qval")
      file <- file[,c("chr","start","end","strand","name","score","fc","output_coverage","pval","qval")]
      file <- unique(makeGRangesFromDataFrame(file,keep.extra.columns=T))
    }
  }
  
  seqlevels(file,pruning.mode="coarse") <- c(paste0("chr",1:22),"chrX","chrY")
  strand(file) <- "*"
  if (remove_mcols) mcols(file) <- NULL
  file
}

paths_to_gr <- function(paths,remove_mcols=FALSE){
  # read in a list of characters representing paths, output single GRanges
  gr_list <- lapply(paths,read_file,remove_mcols=remove_mcols)
  do.call(c,gr_list)
}

lift2hg38 <- function(gr){
  ch <- import.chain("/maps/projects/ralab/people/pcr980/Resource_Genome_annot/hg19ToHg38.over.chain")
  unlist(liftOver(gr, ch))
}

write_narrowPeak <- function(gr,out_path){
  df <- as.data.frame(gr)
  df$width <- NULL
  df$name <- NA
  df$signalValue <- 0
  df$pValue <- 0
  df$qValue <- 0
  df$peak <- 0
  df$score <- 0
  df <- df[,c("seqnames","start","end","name","score","strand","signalValue","pValue","qValue","peak")]
  write.table(df, out_path, quote=F, sep="\t", row.names=F, col.names=F)
}

full_path_to_file_name <- function(path){
  path_segmented <- unlist(strsplit(path,"/"))
  path_segmented[length(path_segmented)]
}

get_genome_info <- function(out_path){
  seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  genomeInfo <- keepStandardChromosomes(seqInfo,species="Homo_sapiens")

  df <- as.data.frame(genomeInfo)
  df$isCircular <- NULL
  df$genome <- NULL
  colnames(df) <- c("end")
  df$end <- df$end -1
  df$start <- 0
  df <- df[,c("start","end")]
  write.table(df,out_path, quote=F, sep="\t", col.names=F)
}



# calculate file name for multilabel file
# single-labeled file name starts with modality
# multi-labeled file name starts with cell type
calculate_multilabel_file_name <- function(name_list){
  # name_list contain 2 file names

  name1_split <- unlist(stringr::str_split(name_list[1],"_"))
  modality1 <- name1_split[1]
  cell1 <- name1_split[2]
  
  name2_split <- unlist(stringr::str_split(name_list[2],"_"))
  modality2 <- name2_split[1]
  cell2 <- name2_split[2]
  
  if (cell1==cell2){
    return(paste(cell1,modality1,modality2,sep="_")) # hepg2_cage_dhs
  } else if (modality1==modality2){
    return(paste(cell1,cell2,modality1,sep="_")) # e.g. hepg2_k562_cage
  } else{
    stop("Invalid input names!")
  }
}


#----------------------------------------------------------
# calculate distances (center to center) and overlaps
#----------------------------------------------------------
calculate_distance <- function(gr1,gr2){
  # report distance from each gr1 to nearest gr2
  # return a vector of distances
  gr1 <- resize(gr1, width=1, fix="center")
  gr2 <- resize(gr2, width=1, fix="center")
  dist <- as.data.frame(distanceToNearest(gr1,gr2,ignore.strand=TRUE))
  dist$distance
}

calculate_distance_self <- function(gr){
  # report distance from each region in gr to nearest region within gr
  # return a vector of distances
  gr <- resize(gr, width=1, fix="center")
  dist <- as.data.frame(distanceToNearest(gr,ignore.strand=TRUE))
  dist$distance
}


full_overlap_info <- function(gr1,gr2){
  # for all possible overlap (represented by index)
  # report length of each overlapping part, width and proportion the of overlapped region
  # return a data frame
  ov_info <- as.data.frame(findOverlaps(gr1,gr2))
  gr1_redun <- gr1[ov_info$queryHits]
  gr2_redun <- gr2[ov_info$subjectHits]
  intersection <- pintersect(gr1_redun,gr2_redun)
  ov_info$gr1_length <- width(gr1_redun)
  ov_info$gr2_length <- width(gr2_redun)
  ov_info$intersection_length <- width(intersection)
  ov_info$ratio <- ov_info$intersection_length/pmin(ov_info$gr1_length,ov_info$gr2_length)
  ov_info
}


calculate_overlap <- function(gr1,gr2){
  # for each region in gr1, what percentage does it overlap with most overlapping region in gr2
  # output maximum of the following two
  # 1. length(overlap)/length(gr1)
  # 2. length(overlap)/length(gr2)
  # return a vector of proportion
  gr1$ratio <- 0
  ov_info <- full_overlap_info(gr1,gr2)
  ov_info <- ov_info %>% 
    group_by(queryHits)%>% 
    filter(ratio == max(ratio)) %>%
    as.data.frame()
  gr1[ov_info$queryHits]$ratio <- ov_info$ratio
  gr1$ratio
}

calculate_overlap_self <- function(gr){
  # for each nearest neighbor pair in gr, what is their percentage of overlap
  # output length(overlap)/min(gr)
  gr <- sortSeqlevels(gr)
  gr <- sort(gr,ignore.strand=TRUE)
  num <- length(gr)
  width(pintersect(gr[1:num-1],gr[2:num]))/pmin(width(gr[1:num-1]),width(gr[2:num]))
}
