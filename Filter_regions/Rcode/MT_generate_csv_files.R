library(rtracklayer)
library(purrr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyr)
source("/isdata/alab/people/pcr980/R/Function_bigwig.R")
source("/isdata/alab/people/pcr980/R/Function_utils.R")
source("/isdata/alab/people/pcr980/R/Function_multilabel.R")
source("/isdata/alab/people/pcr980/R/Function_generate_data_for_ML.R")



# To be very confident:
# POS_OV_THRESH=0.8
# NEG_OV_THRESH=0


#-----------------------------
# specify data path
#----------------------------

dir_base <- "/isdata/alab/people/pcr980/Raw_data/"
universe <- file.path(dir_base,"Bed_processed/universe.bed")

pos_list_k562 <- c(file.path(dir_base,"Bed_processed/CAGE_K562.bed"),
                   file.path(dir_base,"Bed_processed/DHS_K562.bed"),
                   file.path(dir_base,"Bed_processed/STARR_K562.bed"),
                   file.path(dir_base,"Bed_processed/SuRE_K562.bed"))
loose_list_k562 <- c(file.path(dir_base,"Bed_processed/CAGE_K562_loose.bed"),
                     file.path(dir_base,"Bed_processed/DHS_K562_loose.bed"),
                     file.path(dir_base,"Bed_processed/STARR_K562_loose.bed"),
                     file.path(dir_base,"Bed_processed/SuRE_K562_loose.bed"))


pos_list_hepg2 <- c(file.path(dir_base,"Bed_processed/CAGE_HepG2.bed"),
                   file.path(dir_base,"Bed_processed/DHS_HepG2.bed"),
                   file.path(dir_base,"Bed_processed/STARR_HepG2.bed"),
                   file.path(dir_base,"Bed_processed/SuRE_HepG2.bed"))
loose_list_hepg2 <- c(file.path(dir_base,"Bed_processed/CAGE_HepG2_loose.bed"),
                      file.path(dir_base,"Bed_processed/DHS_HepG2_loose.bed"),
                      file.path(dir_base,"Bed_processed/STARR_HepG2_loose.bed"),
                      file.path(dir_base,"Bed_processed/SuRE_HepG2_loose.bed"))



out_dir <- "/isdata/alab/people/pcr980/Dataset_more_distinguishing/"
quantification_type <- "avg"
aug <- T
split_type <- "final_testing"

#----------------------------------------
# Step 1: (option1) generate single-labeled file，together with 8 quantified signal
#----------------------------------------

# HepG2
data_names <- c("cage_hepg2","dhs_hepg2","starr_hepg2","sure_hepg2")
lapply(1:4,function(i){
  bed2csv(pos_list_hepg2[i],
          loose_list_hepg2[i],
          600,
          out_dir,
          data_names[i],
          quantification_type,
          aug,
          split_type)
})

# K562
data_names <- c("cage_k562","dhs_k562","starr_k562","sure_k562")
lapply(1:4,function(i){
  bed2csv(pos_list_k562[i],
          loose_list_k562[i],
          600,
          out_dir,
          data_names[i],
          quantification_type,
          aug,
          split_type)
})
  



#---------------------------
# Merge single labeled file files
#--------------------------


files <- list.files(out_dir,pattern="_train.csv",full.names = T)
merge_files(files,paste0(out_dir,"dat_train.csv"))

files <- list.files(out_dir,pattern="_val.csv",full.names = T)
merge_files(files,paste0(out_dir,"dat_val.csv"))

files <- list.files(out_dir,pattern="_test.csv",full.names = T)
merge_files(files,paste0(out_dir,"dat_test.csv"))


#----------------------------------------
# Step 3: derive multilabel
#----------------------------------------

derive_multilabels(c("cage","dhs","starr","sure"),
                   pos_list_hepg2,loose_list_hepg2,
                   pos_list_k562,loose_list_k562,
                   out_dir,
                   quantification_type,
                   split_type,
                   EDA_path=file.path(dir_base,"EDAs/multilabel_statistics.txt"))



merge_files(c(file.path(out_dir,"hepg2_k562_cage_train.csv"),
              file.path(out_dir,"hepg2_k562_dhs_train.csv"),
              file.path(out_dir,"hepg2_k562_starr_train.csv"),
              file.path(out_dir,"hepg2_k562_sure_train.csv")),
            paste0(out_dir,"dat_distinguishing_train.csv"))

merge_files(c(file.path(out_dir,"hepg2_k562_cage_val.csv"),
              file.path(out_dir,"hepg2_k562_dhs_val.csv"),
              file.path(out_dir,"hepg2_k562_starr_val.csv"),
              file.path(out_dir,"hepg2_k562_sure_val.csv")),
            paste0(out_dir,"dat_distinguishing_val.csv"))











#-----------------------------------------------
# Optional: Add quantified loose positive regions
#-----------------------------------------------
# remove regions already in qualified positive

# HepG2
data_names <- c("cage_hepg2","dhs_hepg2","starr_hepg2","sure_hepg2")
lapply(1:4,function(i){
  bed2csv_loose(pos_list_hepg2[i],
                loose_list_hepg2[i],
                600,
                file.path(dir_base,"CSV_unaugmented",paste0(data_names[i],"_loose.csv"))
  )
})

# K562
data_names <- c("cage_k562","dhs_k562","starr_k562","sure_k562")
lapply(1:4,function(i){
  bed2csv_loose(pos_list_k562[i],
                loose_list_k562[i],
                600,
                file.path(dir_base,"CSV_unaugmented",paste0(data_names[i],"_loose.csv"))
  )
})





#----------------------------------------------------------------------------------
# Step 1: (option 2) generate single-labeled file，together with 1 quantified signal
#----------------------------------------------------------------------------------
cage_hepg2 <- file.path(dir_base,"Bed_processed/CAGE_HepG2.bed")
cage_k562 <- file.path(dir_base,"Bed_processed/CAGE_K562.bed")

quantify_one_rdata_file_for_cage <- function(gr_path,rdata_file,column_name,out_prefix){
  load(rdata_file)
  seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  genomeInfo <- keepStandardChromosomes(seqInfo,species="Homo_sapiens")
  gr <- read_file(gr_path,remove_mcols=T)
  gr$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resize(gr, width=600, fix="center"))
  
  names(gr) <- paste0("Seq",1:length(gr))
  seqinfo(gr) <- genomeInfo
  gr$thick <- resize(gr, width=1, fix="center")
  strand(gr) <- "*"
  expr <- quantifyOlapClusters(CTSSs,clusters=gr)
  expr <- calcTPM(expr)
  df_res <- as.data.frame(log1p(rowMeans(assay(expr,"TPM"))))
  colnames(df_res) <- column_name
  df_res <- df_res[names(gr),,drop=F]
  df_res <- as.data.frame(cbind(as.data.frame(gr),df_res))
  df_res <- df_res[complete.cases(df_res),]
  rownames(df_res) <- NULL
  train_test_split(df_res,"/isdata/alab/people/pcr980/Dataset_CAGE/",out_prefix,"model_selection")
}

cage_hepg2_rdata <- c("cage_hepg2_all_00","cage_hepg2_all_10",
                      "cage_hepg2_C_00","cage_hepg2_C_10",
                      "cage_hepg2_N_00","cage_hepg2_N_10")
lapply(1:6,function(i){
  quantify_one_rdata_file_for_cage(cage_hepg2,
                                   paste0(dir_base,"CAGE_Rdata/",cage_hepg2_rdata[i],".Rdata"),
                                   "cage_hepg2_intensity",
                                   cage_hepg2_rdata[i])
})


cage_k562_rdata <- c("cage_k562_all_00","cage_k562_all_10",
                     "cage_k562_C_00","cage_k562_C_10",
                     "cage_k562_N_00","cage_k562_N_10")
lapply(1:6,function(i){
  quantify_one_rdata_file_for_cage(cage_k562,
                                   paste0(dir_base,"CAGE_Rdata/",cage_k562_rdata[i],".Rdata"),
                                   "cage_k562_intensity",
                                   cage_k562_rdata[i])
})


