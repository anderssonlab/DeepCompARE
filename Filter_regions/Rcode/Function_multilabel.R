library(rtracklayer)
library(purrr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ggplot2)

POS_OV_THRESH=0.8
NEG_OV_THRESH=0

universe <- "/isdata/alab/people/pcr980/Raw_data/Bed_processed/universe.bed"
source("/isdata/alab/people/pcr980/R/Function_utils.R")
source("/isdata/alab/people/pcr980/R/Function_generate_data_for_ML.R")

#----------------------
# Define functions and common data
#---------------------

get_variable_name <- function(var) {
  deparse(substitute(var))
}

report_overlaps <- function(gr_list,data_names){
  # calculate the ovelap from gr_list[[1]] to rest of gr_list
  
  # by default, gr_list[[1]] is the reference
  ovratio_list <- lapply(2:length(gr_list),function(i){
    calculate_overlap(gr_list[[1]],gr_list[[i]])
  })
  ovratio_list <- reduce(ovratio_list,cbind)
  ovratio_df <- as.data.frame(ovratio_list)
  colnames(ovratio_df) <- data_names[2:length(data_names)]
  ovratio_df
}

ov2label <- function(df,df_loose,pos_ov_thresh,neg_ov_thresh){
  # number of datasets is irrelevant
  df[df>=pos_ov_thresh] <- 1
  df[df<=neg_ov_thresh] <- 0
  df[df>neg_ov_thresh & df<pos_ov_thresh] <- NA
  df[df_loose>neg_ov_thresh & df==0] <- NA
  df <- df[complete.cases(df),,drop=F]
  df
}

balance_labels <- function(gr){
  gr$labels <- do.call(paste0,mcols(gr))
  label_table <- table(gr$labels)
  # remove "11", only maintain 01 and 10
  gr <- gr[gr$labels!="11"]
  # get number of samples
  num_sample_per_class <- min(unname(label_table))
  # subsample gr if necessary.
  grs <- lapply(names(label_table),function(this_label){
    gr_sub <- gr[gr$labels==this_label]
    gr_sub$labels <- NULL
    idx <- sample(length(gr_sub),min(num_sample_per_class,length(gr_sub)))
    gr_sub[idx]
  })
  do.call(c,grs)
}


report_labels <- function(gr,EDA_path=NULL){
  labels <- do.call(paste0,mcols(gr))
  my_table <- table(labels)
  print(my_table)
  if (!is.null(EDA_path)){
    write(paste(names(my_table),collapse=" "),EDA_path,append=T)
    write(paste(as.character(my_table),collapse=" "),EDA_path,append=T)
  }
}

calculate_labels_ref_pos<- function(data_names,pos_list,loose_list,my_order,pos_ov_thresh,neg_ov_thresh){
  num_data <- length(data_names)
  # reorder input lists
  data_names <- data_names[my_order]
  pos_list <- pos_list[my_order]
  loose_list <- loose_list[my_order]
  # use pos as ref,not loose
  pos_list[[1]] <- resize(pos_list[[1]],width=600,fix="center")
  loose_list[[1]] <- pos_list[[1]]
  
  df_pos <- report_overlaps(pos_list,data_names)
  colnames(df_pos) <- data_names[2:num_data]
  df_loose <- report_overlaps(loose_list,data_names)
  colnames(df_loose) <- data_names[2:num_data]
  rnames <- union(rownames(df_pos),rownames(df_loose))
  df_pos <- df_pos[rnames,,drop=F]
  df_loose <- df_loose[rnames,,drop=F]
  df_pos <- ov2label(df_pos,df_loose,pos_ov_thresh,neg_ov_thresh)
  rnames <- rownames(df_pos)
  
  gr_ref <- pos_list[[1]]
  gr_ref <- gr_ref[as.numeric(rnames)]
  
  mcols(gr_ref)[[data_names[1]]] <- 1
  mcols(gr_ref)[[data_names[2]]] <- df_pos[,1]
  
  if (num_data>2) mcols(gr_ref)[[data_names[3]]] <- df_pos[,2]
  if (num_data>3) mcols(gr_ref)[[data_names[4]]] <- df_pos[,3]
  
  gr_ref
}


calculate_labels_pos <- function(data_names,pos_list,loose_list,pos_ov_thresh,neg_ov_thresh){
  # returns gr with at least one positive label. Each file in pos_list gets to be the reference once.
  if (length(data_names)==4) order_list <- list(c(1,2,3,4),c(2,1,3,4),c(3,1,2,4),c(4,1,2,3))
  if (length(data_names)==3) order_list <- list(c(1,2,3),c(2,1,3),c(3,1,2))
  if (length(data_names)==2) order_list <- list(c(1,2),c(2,1))
  gr <- lapply(order_list,function(this_order){
    calculate_labels_ref_pos(data_names,pos_list,loose_list,this_order,pos_ov_thresh,neg_ov_thresh)
  })
  gr <- do.call(c,gr)
  unique(gr)
}



derive_labels <- function(data_names,pos_list,loose_list,out_dir,out_prefix,
                          quantification_type,split_type,EDA_path=NULL,
                          pos_ov_thresh=POS_OV_THRESH,neg_ov_thresh=NEG_OV_THRESH){
  # derive all possible labels for given pos_list
  # data_names may have length 2,3,4
  
  # read in all BED files and store in GRangesList object
  pos_list <- lapply(pos_list,read_file,remove_mcols = T)
  loose_list <- lapply(loose_list,read_file,remove_mcols = T)

  # use positive as reference 
  gr <- calculate_labels_pos(data_names,pos_list,loose_list,pos_ov_thresh,neg_ov_thresh)
  
  # balance labels
  if (!is.null(EDA_path)) write(paste0(out_prefix," before balancing label"),EDA_path,append=T)
  report_labels(gr,EDA_path)
  gr <- balance_labels(gr)
  if (!is.null(EDA_path)) write(paste0(out_prefix," after balancing label"),EDA_path,append=T)
  report_labels(gr,EDA_path)
  
  # RC augmentation
  strand(gr) <- "+"
  gr_aug <- gr
  strand(gr_aug) <- "-"
  gr <- c(gr,gr_aug)
  
  # get sequence
  message("Getting sequence of 600 bp")
  gr$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr)
  
  message("Adding signal intensity")
  names(gr) <- paste0("Seq",1:length(gr))
  
  if (quantification_type=="avg"){
    df_quantification <- quantify_signals(resize(gr, width=600, fix="center"),"/isdata/alab/people/pcr980/Raw_data/BW_processed")
  } else if (quantification_type=="by_nucleosome"){
    df_quantification <- quantify_3_nucleosomes(resize(gr, width=600, fix="center"),"/isdata/alab/people/pcr980/Raw_data/BW_processed")
  } else {
    stop("quantifycation type not recoginzed")
  }
  
  df_quantification <- log1p(df_quantification)
  df <- as.data.frame(cbind(as.data.frame(gr),df_quantification))
  df <- df[complete.cases(df),]
  # write output file
  message("Performing train-test spliting")
  train_test_split(df,out_dir,out_prefix,split_type)
}












derive_multilabels <- function(data_names,
                               pos_list_hepg2,loose_list_hepg2,
                               pos_list_k562,loose_list_k562,
                               out_dir,
                               quantification_type,
                               split_type,
                               EDA_path){

  # step 4: Derive 2 labels of 1 modality of 2 cell types
  lapply(1:length(data_names),function(i){
    write(paste0("hepg2_k562_",data_names[i]),EDA_path,append=T)
    derive_labels(c(paste0(data_names[i],"_hepg2_class"),paste0(data_names[i],"_k562_class")),
                  c(pos_list_hepg2[i],pos_list_k562[i]),
                  c(loose_list_hepg2[i],loose_list_k562[i]),
                  out_dir,
                  paste0("hepg2_k562_",data_names[i]),
                  quantification_type,
                  split_type,
                  EDA_path
    )
  })
}












#---------------
# Trash bin
#---------------


calculate_labels_neg<- function(universe,data_names,loose_list){
  # returns gr with label 0000/000/00
  universe <- read_file(universe,remove_mcols = T)
  universe <- resize(universe, width=600, fix="center")
  gr_loose <- do.call(c,loose_list)
  ov <- as.data.frame(findOverlaps(universe,gr_loose,ignore.strand=TRUE))
  neg <- universe[-unique(ov$queryHits)]
  for (name in data_names) mcols(neg)[[name]] <- 0
  neg
}
