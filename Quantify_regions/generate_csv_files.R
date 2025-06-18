setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions")
library(rtracklayer)
library(purrr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyr)
source("Function_bigwig.R")
source("Function_utils.R")
source("Function_multilabel.R")
source("Function_generate_data_for_ML.R")



# To be very confident:
# POS_OV_THRESH=0.8
# NEG_OV_THRESH=0


#-----------------------------
# specify data path
#----------------------------

dir_bed <- "/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/"
universe <- file.path(dir_bed,"universe.bed")

pos_list_k562 <- c(file.path(dir_bed,"CAGE_K562.bed"),
                   file.path(dir_bed,"DHS_K562.bed"),
                   file.path(dir_bed,"STARR_K562.bed"),
                   file.path(dir_bed,"SuRE_K562.bed"))
loose_list_k562 <- c(file.path(dir_bed,"CAGE_K562_loose.bed"),
                     file.path(dir_bed,"DHS_K562_loose.bed"),
                     file.path(dir_bed,"STARR_K562_loose.bed"),
                     file.path(dir_bed,"SuRE_K562_loose.bed"))


pos_list_hepg2 <- c(file.path(dir_bed,"CAGE_HepG2.bed"),
                   file.path(dir_bed,"DHS_HepG2.bed"),
                   file.path(dir_bed,"STARR_HepG2.bed"),
                   file.path(dir_bed,"SuRE_HepG2.bed"))
loose_list_hepg2 <- c(file.path(dir_bed,"CAGE_HepG2_loose.bed"),
                      file.path(dir_bed,"DHS_HepG2_loose.bed"),
                      file.path(dir_bed,"STARR_HepG2_loose.bed"),
                      file.path(dir_bed,"SuRE_HepG2_loose.bed"))




#-------------------------------------------------------------------------------
# Step 1: generate single-labeled file，together with 8 quantified signal
#-------------------------------------------------------------------------------

# Split into 5 chunks, HepG2
generate_single_labeled_dataset(pos_list_hepg2,loose_list_hepg2,
                                600,
                                out_dir="CV5_data_chunks/",
                                data_names=c("cage_hepg2","dhs_hepg2","starr_hepg2","sure_hepg2"),
                                quantification_type="avg",
                                aug=T,
                                split_type="cv5")

# Split into 5 chunks, K562
generate_single_labeled_dataset(pos_list_k562,loose_list_k562,
                                600,
                                out_dir="CV5_data_chunks/",
                                data_names=c("cage_k562","dhs_k562","starr_k562","sure_k562"),
                                quantification_type="avg",
                                aug=T,
                                split_type="cv5")


# Final testing, HepG2

generate_single_labeled_dataset(pos_list_hepg2,loose_list_hepg2,
                                600,
                                out_dir="/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/",
                                data_names=c("cage_hepg2","dhs_hepg2","starr_hepg2","sure_hepg2"),
                                quantification_type="avg",
                                aug=T,
                                split_type="train_val_test")

# Final testing, K562
generate_single_labeled_dataset(pos_list_k562,loose_list_k562,
                                600,
                                out_dir="/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/",
                                data_names=c("cage_k562","dhs_k562","starr_k562","sure_k562"),
                                quantification_type="avg",
                                aug=T,
                                split_type="train_val_test")

#----------------------------------------
# Step 2: derive multilabel
#----------------------------------------
# Split into 5 chunks
derive_multilabels(c("cage","dhs","starr","sure"), # RC augmentation by default
                   pos_list_hepg2,loose_list_hepg2,
                   pos_list_k562,loose_list_k562,
                   out_dir="/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions/CV5_data_chunks/",
                   quantification_type="avg",
                   split_type="cv5",
                   EDA_path="EDAs/multilabel_statistics.txt")


# Final testing
derive_multilabels(c("cage","dhs","starr","sure"), # RC augmentation by default
                   pos_list_hepg2,loose_list_hepg2,
                   pos_list_k562,loose_list_k562,
                   out_dir="/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/",
                   quantification_type="avg",
                   split_type="train_val_test",
                   EDA_path="EDAs/multilabel_statistics.txt")



#---------------------------
# Step 3: Merge
#--------------------------
# Merge chunks into cross-validation data
merge_5cv("CV5_data_chunks/",
          "/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_5cv_with_multilabel",
          with_multilabel = T)
merge_5cv("CV5_data_chunks/",
          "/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_5cv_without_multilabel",
          with_multilabel = F)




out_dir <- "/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/"

files <- list.files(out_dir,pattern="_train.csv",full.names = T)
merge_files(files,paste0(out_dir,"dat_train.csv"))

files <- list.files(out_dir,pattern="_val.csv",full.names = T)
merge_files(files,paste0(out_dir,"dat_val.csv"))

files <- list.files(out_dir,pattern="_test.csv",full.names = T)
merge_files(files,paste0(out_dir,"dat_test.csv"))







#-----------------------------------------------
# Step 4: Quantify loose positive regions for pure regression
#-----------------------------------------------
# remove regions already in qualified positive

# HepG2
generate_loose_dataset(
  pos_list_hepg2,
  loose_list_hepg2,
  600,
  "CV5_loose_chunks/",
  data_names=c("cage_hepg2","dhs_hepg2","starr_hepg2","sure_hepg2")
)


# K562
generate_loose_dataset(
  pos_list_k562,
  loose_list_k562,
  600,
  "CV5_loose_chunks/",
  data_names=c("cage_k562","dhs_k562","starr_k562","sure_k562")
)


# merge into 5-fold cross validation datasets
merge_5cv("CV5_loose_chunks/",
          "/maps/projects/ralab/people/pcr980/DeepCompare/Datasets/Dataset_5cv_loose",
          with_multilabel = F)






