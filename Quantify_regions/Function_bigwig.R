source("Function_utils.R")
library(CAGEfightR)
library(rtracklayer)
library(liftOver)
library(purrr)


# add or average bigwig signals.
combine_bw <- function(dir,flist,outname,method,num_samples=NULL){
  cov_list <- lapply(flist,function(fname){
    bw_file <- import.bw(file.path(dir,fname))
    coverage(bw_file,weight=bw_file$score)
  })
  combined_cov <- cov_list[[1]]
  for (item in cov_list[2:length(flist)]){
    combined_cov <- combined_cov+item
  }
  if (method=="avg"){
    if (num_samples=="infer") combined_cov <- combined_cov/length(flist)
    if (is.numeric(num_samples))  combined_cov <- combined_cov/num_samples
  }
  gr <- as(combined_cov,"GRanges")
  export.bw(gr,outname)
}


quantifyOlapClusters <- function (object, clusters, inputAssay = "counts", sparse = FALSE) 
{
  assertthat::assert_that(methods::is(object, "RangedSummarizedExperiment"), 
                          assertthat::not_empty(object), 
                          isDisjoint(object), 
                          methods::is(clusters,"GRanges"), 
                          assertthat::not_empty(clusters), 
                          identical(seqinfo(object), seqinfo(clusters)),
                          inputAssay %in% assayNames(object)
                          )
  message("Finding overlaps...")
  hits <- findOverlaps(query = object, subject = clusters, 
                       select = "all")
  keep <- queryHits(hits)
  cluster_hits <- factor(subjectHits(hits), levels = seq_along(clusters))
  missing_hits <- is.na(hits)
  if (all(missing_hits)) {
    warning("The supplied clusters had no overlapping CTSSs!")
  }
  message("Aggregating within clusters...")
  mat <- utilsAggregateRows(x = assay(object, inputAssay)[keep,,drop=F], 
                            group = cluster_hits, drop = FALSE, sparse = sparse)
  stopifnot(nrow(mat) == length(clusters))
  rownames(mat) <- names(clusters)
  o <- SummarizedExperiment(assays = SimpleList(mat), rowRanges = clusters, 
                            colData = colData(object))
  assayNames(o) <- inputAssay
  o
}



# given GRanges, find average bigwig signal on this range
quantify_one_bw_file <- function(gr,bw_file){
  # gr must be named
  message("Quantify on ",bw_file)
  strand(gr) <- "*"
  bw_score <- import.bw(bw_file,which=gr)
  ov_info <- as.data.frame(findOverlaps(gr,bw_score))
  ov_info$width <- width(bw_score[ov_info$subjectHits])
  ov_info$score <- bw_score$score[ov_info$subjectHits]
  ov_info$score_sum <- ov_info$width*ov_info$score
  ov_info$seq_id <- names(gr)[ov_info$queryHits]
  df_res <- aggregate(ov_info$score_sum,by=list(ov_info$seq_id),sum)
  colnames(df_res) <- c("seq_id",full_path_to_file_name(bw_file))
  rownames(df_res) <- df_res$seq_id
  df_res$seq_id <- NULL
  df_res[,1] <- df_res[,1]/width(gr[rownames(df_res)])
  colnames(df_res) <- gsub(".bw","_intensity",colnames(df_res))
  colnames(df_res) <- gsub(".bigWig","_intensity",colnames(df_res))
  df_res
}


# given GRanges, find average cage TPM expression on this range
quantify_one_rdata_file <- function(gr,rdata_file){
  load(rdata_file)
  message("Quantify on ",rdata_file)
  seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  genomeInfo <- keepStandardChromosomes(seqInfo,species="Homo_sapiens")
  seqinfo(gr) <- genomeInfo
  strand(gr) <- "*"
  expr <- suppressMessages(quantifyOlapClusters(CTSSs,clusters=gr))
  expr <- suppressMessages(calcTPM(expr))
  df_res <- as.data.frame(rowMeans(assay(expr,"TPM")))
  colnames(df_res) <- gsub(".Rdata","_intensity",full_path_to_file_name(rdata_file))
  df_res
}

# given GRanges, find average expression on this range according to 8 bw tracks
quantify_signals <- function(gr,dir_bw){
  # gr must be named

  files_bw <- list.files(dir_bw,full.names = T)
  annotations <- lapply(files_bw,function(fname){
    if (grepl(".bw",fname)){
      annotation <- quantify_one_bw_file(gr,fname)
    } else if (grepl(".bigWig",fname)){
      annotation <- quantify_one_bw_file(gr,fname)
    } else if (grepl(".Rdata",fname)){
      annotation <- quantify_one_rdata_file(gr,fname)
    } else{
      stop("Signal file format unrecognized.")
    }
    annotation <- annotation[names(gr),,drop=F]
    annotation
  })
  df_res <- do.call(cbind,annotations)
  stopifnot(rownames(df_res)==names(gr))
  df_res
}


quantify_3_nucleosomes <- function(gr,dir_bw){
  # gr must be named
  gr_left <- resize(gr,width=200,fix="start")
  gr_middle <- resize(gr,width=200,fix="center")
  gr_right <- resize(gr,width=200,fix="end")
  
  df_left <- quantify_signals(gr_left,dir_bw)
  df_middle <- quantify_signals(gr_middle,dir_bw)
  df_right <- quantify_signals(gr_right,dir_bw)
  
  colnames(df_left) <- paste0(colnames(df_left),"_left")
  colnames(df_middle) <- paste0(colnames(df_middle),"_middle")
  colnames(df_right) <- paste0(colnames(df_right),"_right")
  
  as.data.frame(cbind(df_left,df_middle,df_right))
  
  
}


