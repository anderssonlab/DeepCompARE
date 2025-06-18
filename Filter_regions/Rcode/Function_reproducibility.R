
source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/Generate_data/Function_utils.R")
library(rtracklayer)
library(gridExtra)
library(grid)
library(gridBase)
#-------------------------------
# Reproducibility EDA
#-------------------------------

reproducibility_within_modality_EDA <- function(data_dir,thresh,title_prefix){
  # input data directory, output reproducibility profile on HepG2 and K562 separately
  fhepg2 <- list.files(data_dir,pattern="hepg2",ignore.case=T,full.name=T)
  fk562 <- list.files(data_dir,pattern="k562",ignore.case=T,full.name=T)
  par(mfrow=c(2,1))
  reproducibility_EDA_single(fhepg2,thresh,paste0(title_prefix," HepG2"))
  reproducibility_EDA_single(fk562,thresh,paste0(title_prefix," K562"))
}

reproducibility_EDA_single <- function(fnames,thresh,title){
  gr <- paths_to_gr(fnames,remove_mcols=T)
  gr_cov <- coverage(gr)
  gr <- as(gr_cov,"GRanges")
  gr$length <- width(gr)
  df <- as.data.frame(gr)
  df <- df[df$score>0,]

  df$score <- as.factor(df$score)
  p1 <- ggplot(data=df,aes(x=length,color=score))+
    geom_density()+
    ggtitle(paste0(title," Length v.s. score"))+
    geom_vline(xintercept=thresh)
  
  df$score <- as.numeric(df$score)
  p2 <- hist(df$score,plot = FALSE)
  
  df <- df[df$length>thresh,]
  p3 <- hist(df$score,plot = FALSE)
  
  layout(matrix(1:2, ncol = 2))
  print(p1)
  plot(p2, col = 'lightblue', main=paste0(title," Score distribution"))
  plot(p3, col = 'lightblue', main=paste0(title," Score above length thresh"))
}



#----------------------------------
# reproducibility within modality
#----------------------------------

export_reproducible_single <- function(bed_files,thresh_length,thresh_score,out_name,lift){
  bed_regions <- paths_to_gr(bed_files,remove_mcols=T)
  gr_rep <- as(coverage(bed_regions),"GRanges")
  gr_rep$length <- width(gr_rep)
  gr_rep <- gr_rep[gr_rep$score>=thresh_score & gr_rep$length>=thresh_length]
  if (lift) gr_rep <- lift2hg38(gr_rep)
  message("Output to ",out_name)
  export.bed(gr_rep,out_name)
}


export_reproducible_within_modality <- function(bed_dir,thresh_length,thresh_score,out_prefix,lift=F){
  # input data directory, output reproducible regions on HepG2 and K562 separately
  
  bed_files_hepg2 <- list.files(bed_dir,pattern="hepg2",ignore.case=T,full.names=T)
  export_reproducible_single(bed_files_hepg2,
                             thresh_length,
                             thresh_score,
                             paste0(out_prefix,"_HepG2.bed"),
                             lift
  )
  
  bed_files_k562 <- list.files(bed_dir,pattern="k562",ignore.case=T,full.names=T)
  export_reproducible_single(bed_files_k562,
                             thresh_length,
                             thresh_score,
                             paste0(out_prefix,"_K562.bed"),
                             lift
  )
}

export_loose_within_modality <- function(dir,prefix){
  fhepg2 <- list.files(dir,pattern="hepg2",ignore.case = T,full.names = T)
  fk562 <- list.files(dir,pattern="k562",ignore.case = T,full.names = T)
  loose_hepg2 <- paths_to_gr(fhepg2)
  loose_k562 <- paths_to_gr(fk562)
  export.bed(loose_hepg2,paste0("/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/",prefix,"_HepG2_loose.bed"))
  export.bed(loose_k562,paste0("/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/",prefix,"_K562_loose.bed"))
}


#-------------------------------------
# Reproducibility between modalities
#-------------------------------------
meaningfully_intersected <- function(gr1,gr2,intersect_thresh){
  # output subset of gr1 that includes at least "intersect_thresh" of either gr1 or gr2
  
  ov_info <- full_overlap_info(gr1,gr2)
  ov_info <- ov_info[ov_info$ratio>=intersect_thresh,]
  idx <- unique(ov_info$queryHits)
  idx <- idx[!is.na(idx)]
  gr1[idx]
}


single_file_supported_across_modality <- function(file_query, file_reproducible_within_modality, files_rest){
  # given a single file
  # output regions that are not supported within modality but supported across modalities
  
  # read in query file
  gr_query <- read_file(file_query,remove_mcols=T)
  # remove regions reproducible within modality
  gr_reproducible_within_modality <- read_file(file_reproducible_within_modality,remove_mcols=T)
  ov <- findOverlaps(gr_query,gr_reproducible_within_modality)
  gr_query <- gr_query[-unique(queryHits(ov))]
  # search for support across modality
  gr_reproducible_across_modalities <- lapply(files_rest,function(file){
    gr2 <- read_file(file,remove_mcols=T)
    meaningfully_intersected(gr_query,gr2,0.8)
  })
  gr_reproducible_across_modalities <- do.call(c,gr_reproducible_across_modalities)
  unique(gr_reproducible_across_modalities)
}

single_celltype_supported_across_modality <- function(dir_irreproducible,dir_reproducible,cell){
  # given a cell type, a specific modality,
  # output regions that are not supported within modality but supported across modalities
  
  files_query <- list.files(dir_irreproducible,pattern=cell,ignore.case = T,full.names = T)
  
  file_reproducible_within_modality <- list.files(dir_reproducible,pattern=cell,ignore.case=T,full.names = T)
  
  dirs_rest <- dirs_irreproducible[dirs_irreproducible!=dir_irreproducible]
  files_rest <- lapply(dirs_rest,function(dir) list.files(dir,pattern=cell,ignore.case=T,full.names = T))
  files_rest <- unlist(files_rest)
  
  gr_rep <-lapply(files_query,function(file_query){
    gr <- single_file_supported_across_modality(file_query, file_reproducible_within_modality, files_rest)
    message(paste0(file_query,":",length(gr)," regions rescued."))
    gr
  })
  gr_rep <- do.call(c,gr_rep)
  gr_rep <- unique(gr_rep)
}


combine_peaks_reproducible_and_rescued <- function(dir,cell,out_name){
  f1 <- list.files(paste0(dir,"Peaks_reproducible_within_modality/"),pattern=cell,ignore.case=T,full.names = T)
  f2 <- list.files(paste0(dir,"Peaks_rescued/"),pattern=cell,ignore.case=T,full.names = T)
  gr1 <- read_file(f1,remove_mcols=T)
  gr2 <- read_file(f2,remove_mcols=T)
  gr <- unique(c(gr1,gr2))
  export.bed(gr,out_name)
}