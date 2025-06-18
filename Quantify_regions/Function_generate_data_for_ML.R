library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
source("Function_utils.R")
source("Function_bigwig.R")


bed2csv <- function(pos_file,pos_loose_file,my_width,out_dir,out_prefix,quantification_type,aug,split_type){
  # Write sequences together with binary label and signal value for one file (e.g. CAGE_HepG2)
  
  message("Processing",pos_file)
  
  # get universe
  universe <-read_file("/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/universe.bed",remove_mcols=T)
  
  # get negatives by removing loose positives
  pos_loose <- read_file(pos_loose_file,remove_mcols=T) 
  ov <- as.data.frame(findOverlaps(universe,pos_loose,ignore.strand=TRUE))
  neg <- universe[-unique(ov$queryHits)]
  
  # get positives
  pos <- read_file(pos_file,remove_mcols=T) 
  
  # remove metadata, add strand information
  mcols(neg) <- NULL
  mcols(pos) <- NULL
  strand(pos) <- "+"
  strand(neg) <- "+"
  
  # downsample neg
  message(length(pos), " positive samples found. Downsample negatives accordingly.")
  if (length(neg)>length(pos)){
    idx <- sample(length(neg),length(pos))
    neg <- neg[idx]
  }
  
  # include binary label
  colname <- tolower(full_path_to_file_name(pos_file))
  colname <- gsub(".bed","",colname)
  message("Adding binary label at column ",paste0(colname,"_class"))
  mcols(pos)[[paste0(colname,"_class")]] <- 1
  mcols(neg)[[paste0(colname,"_class")]] <- 0
  
  # merge positive and negative
  gr <- c(pos,neg)

  # get sequence (always 600bp)
  message("Getting sequence of 600 bp")
  gr$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resize(gr, width=my_width, fix="center"))
  
  # data augmentation if necessary
  if (aug){
    message("Performing data augmentation.")
    gr_aug <- gr
    strand(gr_aug) <- "-"
    gr_aug$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resize(gr_aug, width=my_width, fix="center"))
    gr <- c(gr,gr_aug)
    rm(gr_aug)
  }
  
  # include quantified signal intensity
  message("Adding signal intensity")
  names(gr) <- paste0("Seq",1:length(gr))
  if (quantification_type=="avg"){
    df_quantification <- quantify_signals(resize(gr, width=my_width, fix="center"),
                                          "/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions/Quant_signal_processed")
  } else if (quantification_type=="by_nucleosome"){
    df_quantification <- quantify_3_nucleosomes(resize(gr, width=my_width, fix="center"),
                                                "/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions/Quant_signal_processed")
  } else {
    stop("quantification type not recoginzed")
  }
  
  df_quantification <- log1p(df_quantification)
  
  df <- as.data.frame(cbind(as.data.frame(gr),df_quantification))
  df <- df[complete.cases(df),]
  rownames(df) <- NULL
  message("Performing data spliting") 
  message("Write to ",out_dir)
  
  if (split_type=="cv5"){
    cv5_split(df,out_dir,out_prefix)
  } else if (split_type=="train_val_test"){
    train_test_split(df,out_dir,out_prefix)
  } else {
    stop("split type", split_type, "not recognized.")
  }
}


generate_single_labeled_dataset <- function(pos_list,
                                            loose_list,
                                            my_width,
                                            out_dir,
                                            data_names,
                                            quantification_type,
                                            aug,
                                            split_type){
  # a wrap-up function to apply bed2csc to every file in the file list
  lapply(1:length(pos_list),function(i){
    bed2csv(pos_list[i],
            loose_list[i],
            my_width,
            out_dir,
            data_names[i],
            quantification_type,
            aug,
            split_type)
  })
}




cv5_split <- function(df,out_dir,out_prefix){
  # remove chr2/chr3 for final testing
  df <- df[!df$seqnames %in% c("chr2","chr3"),]
  # split the rest data into 5 chunks
  n_samples <- dim(df)[1]
  idx <- sample(n_samples,n_samples)
  break1 <- ceiling(n_samples*0.2)
  break2 <- ceiling(n_samples*0.4)
  break3 <- ceiling(n_samples*0.6)
  break4 <- ceiling(n_samples*0.8)
  write.csv(df[1:break1,],paste0(out_dir,out_prefix,"_1.csv"),row.names = F)
  write.csv(df[break1:break2,],paste0(out_dir,out_prefix,"_2.csv"),row.names = F)
  write.csv(df[break2:break3,],paste0(out_dir,out_prefix,"_3.csv"),row.names = F)
  write.csv(df[break3:break4,],paste0(out_dir,out_prefix,"_4.csv"),row.names = F)
  write.csv(df[break4:n_samples,],paste0(out_dir,out_prefix,"_5.csv"),row.names = F)
}




train_test_split <- function(df,out_dir,out_prefix){
    # get test data
    idx <- which(df$seqnames %in% c("chr2","chr3"))
    df_test <- df[idx,]
    df <- df[-idx,]
    # train-validation split
    n_samples <- dim(df)[1]
    idx <- sample(n_samples,n_samples*0.02)
    df_train <- df[-idx,]
    df_val <- df[idx,]

    message("Write file ",out_dir,out_prefix,"_train.csv")
    write.csv(df_train,paste0(out_dir,out_prefix,"_train.csv"),row.names = F)
    message("Write file ",out_dir,out_prefix,"_val.csv")
    write.csv(df_val,paste0(out_dir,out_prefix,"_val.csv"),row.names = F)
    message("Write file ",out_dir,out_prefix,"_test.csv")
    write.csv(df_test,paste0(out_dir,out_prefix,"_test.csv"),row.names = F)
}

  
merge_files <- function(files,out_name){
  df_lists <- lapply(files, function(f) {
    read.csv(f)
  })
  df <- reduce(df_lists,rbind.fill)
  
  df[is.na(df)] <- -1
  dup_idx <- which(duplicated(df[,"seq"]))
  if (length(dup_idx)>0){
    df <- df[-dup_idx,]
  }
  write.csv(df,out_name)
}



merge_5cv <- function(dir_data_chunk,out_dir,with_multilabel){
  # merge 5 data chunks into 5 training-validation datasets for 5-fold cross validation
  all_files <- list.files(dir_data_chunk,full.names = T)
  
  if (!with_multilabel){
    all_files <- all_files[!grepl("hepg2_k562",all_files)]
  }
  
  lapply(c("_1","_2","_3","_4","_5"), function(suffix){ 
    files_train <- all_files[!grepl(suffix,all_files)]
    files_val <- all_files[grepl(suffix,all_files)]
    dir.create(paste0(out_dir,suffix))
    merge_files(files_train,paste0(out_dir,suffix,"/dat_train.csv"))
    merge_files(files_val,paste0(out_dir,suffix,"/dat_val.csv"))
  })
}





bed2csv_loose <- function(pos_file,pos_loose_file,my_width,out_dir,out_prefix,aug=T,split_type="cv5"){
  # quantify union of loose and positive regions
  # don't add binary labels
  # for regression task alone.
  message("Processing",pos_loose_file)
  
  # get positives
  pos <- read_file(pos_file,remove_mcols=T) 
  
  # combine loose positives with positive
  # if overlap, select positive first.
  loose <- read_file(pos_loose_file,remove_mcols=T) 
  ov <- as.data.frame(findOverlaps(loose,pos,ignore.strand=TRUE))
  loose <- loose[-unique(ov$queryHits)]
  loose <- c(pos,loose)
  loose <- resize(loose,width=my_width,fix="center")
  
  # remove regions overlapping blacklist for SuRE-Seq
  if (grepl("sure",out_prefix)){
    message("Remove regions overlapping SNP")
    black_list <- read_file("/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/SuRE_SNP_hg38.bed",remove_mcols=T) 
    ov <- as.data.frame(findOverlaps(loose,black_list,ignore.strand=TRUE))
    loose <- loose[-unique(ov$queryHits)]
  }
  
  # remove metadata
  mcols(loose) <- NULL
  
  # get sequence (always 600bp)
  message("Getting sequence of 600 bp")
  strand(loose) <- "+"
  loose$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, loose)
  
  # data augmentation if necessary
  if (aug){
    message("Performing data augmentation.")
    loose_aug <- loose
    strand(loose_aug) <- "-"
    loose_aug$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resize(loose_aug, width=my_width, fix="center"))
    loose <- c(loose,loose_aug)
  }
  
  message(length(loose), " positive samples found.")
  names(loose) <- paste0("Seq",1:length(loose))
  
  # include quantified signal intensity
  message("Adding signal intensity")
  df_quantification <- quantify_signals(resize(loose, width=my_width, fix="center"),"/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions/Quant_signal_processed/")
  df_quantification <- log1p(df_quantification)
  df <- as.data.frame(cbind(as.data.frame(loose),df_quantification))
  df <- df[complete.cases(df),]
  
  message("Write to ",out_dir)
  
  if (split_type=="cv5"){
    cv5_split(df,out_dir,out_prefix)
  } else if (split_type=="train_val_test"){
    train_test_split(df,out_dir,out_prefix)
  } else {
    stop("split type", split_type, "not recognized.")
  }
}



generate_loose_dataset <- function(pos_list,loose_list,my_width,out_dir,data_names){
  # quantify union of loose and positive regions
  # don't add binary labels
  # for regression task alone.
  lapply(1:4,function(i){
    bed2csv_loose(pos_list[i],
                  loose_list[i],
                  600,
                  out_dir,
                  data_names[i])
  })
}



