source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/Generate_data/Function_utils.R")
library(patchwork)



rdhs <- unique(import("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/rDHS_ENCFF619EJE.bed",format="BED"))


single_file_validity_EDA <- function(f){
  # plot region width, dist2nearest within file, dist2nearest rDHS(ENCFF619EJE.bed), intersection/union rDHS
  par(mfrow=c(1,4))
  gr <- read_file(f,remove_mcols = T)
  fname <- full_path_to_file_name(f)
  fname <- gsub(".bed","",fname)
  fname <- gsub(".narrowPeak","",fname)
  hist(width(gr),
       main=paste0("Width ", fname),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(log2(calculate_distance_self(gr)),
       main=paste0("Dist2Nearest ", fname),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(log2(calculate_distance(gr,rdhs)),
       main=paste0("Dist2rdhs ", fname),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  rdhs <- subsetByOverlaps(rdhs,gr)
  hist(calculate_overlap(gr,rdhs),
       main=paste0(fname," with intersecting rDHS"),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
}

single_modality_validity_EDA <- function(dir){
  # plot overlaps/distance between cell types within modality
  message("Process ", dir)
  files <- list.files(dir,full.names = T)
  gr1 <- read_file(files[1],remove_mcols = T) 
  gr2 <- read_file(files[2],remove_mcols = T)
  fname1 <- full_path_to_file_name(files[1])
  fname2 <- full_path_to_file_name(files[2])
  par(mfrow=c(1,4))
  hist(calculate_overlap(gr1,gr2),
       main=paste0("intersection ",fname1," ",fname2),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(calculate_overlap(gr2,gr1),
       main=paste0("intersection ",fname2," ",fname1),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(log2(calculate_distance(gr1,gr2)),
       main=paste0("Dist from ", fname1, " to nearest ", fname2),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  hist(log2(calculate_distance(gr2,gr1)),
       main=paste0("Dist from ", fname2, " to nearest ", fname1),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  
}

modality_compatibility_single_cell_type_EDA <- function(files){
  # plot overlaps/distance between modalities within cell type
  file_order <- list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  lapply(1:6,function(i){
    this_order <- file_order[[i]]
    f1 <- files[this_order[1]]
    f2 <- files[this_order[2]]
    gr1 <- read_file(f1,remove_mcols = T) 
    gr2 <- read_file(f2,remove_mcols = T) 
    f1 <- full_path_to_file_name(f1)
    f2 <- full_path_to_file_name(f2)
    par(mfrow=c(1,4))
    hist(calculate_overlap(gr1,gr2),
         main=paste0("intersection ",f1," ",f2),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    hist(calculate_overlap(gr2,gr1),
         main=paste0("intersection ",f2," ",f1),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    hist(log2(calculate_distance(gr1,gr2)),
         main=paste0("Dist from ", f1, " to nearest ", f2),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    hist(log2(calculate_distance(gr2,gr1)),
         main=paste0("Dist from ", f2, " to nearest ", f1),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  })
  dev.off()
}


modality_compatibility_single_celltype_EDA <- function(files){
  # plot overlaps/distance between modalities within cell type
  file_order <- list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  lapply(1:6,function(i){
    this_order <- file_order[[i]]
    f1 <- files[this_order[1]]
    f2 <- files[this_order[2]]
    gr1 <- read_file(f1,remove_mcols = T) 
    gr2 <- read_file(f2,remove_mcols = T) 
    f1 <- full_path_to_file_name(f1)
    f2 <- full_path_to_file_name(f2)
    par(mfrow=c(1,4))
    hist(calculate_overlap(gr1,gr2),
         main=paste0("intersection ",f1," ",f2),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    hist(calculate_overlap(gr2,gr1),
         main=paste0("intersection ",f2," ",f1),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    hist(log2(calculate_distance(gr1,gr2)),
         main=paste0("Dist from ", f1, " to nearest ", f2),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    hist(log2(calculate_distance(gr2,gr1)),
         main=paste0("Dist from ", f2, " to nearest ", f1),
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  })
}

modality_compatibility_EDA <- function(dir){
  fhepg2 <- list.files(dir,pattern="hepg2",ignore.case=T,full.names=T)
  fk562 <- list.files(dir,pattern="k562",ignore.case=T,full.names=T)
  modality_compatibility_single_celltype_EDA(fhepg2)
  modality_compatibility_single_celltype_EDA(fk562)
}