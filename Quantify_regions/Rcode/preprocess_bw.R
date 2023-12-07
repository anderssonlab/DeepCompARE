setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Quantify_regions")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Scripts_R/Generate_data/Function_bigwig.R")

#------------------------------------
# average bigwig files of replicates
#-----------------------------------
out_dir <- "Quant_signal_processed/"


process_bw_wrapper <- function(dir,out_prefix,num_samples="infer"){
  flist_hepg2 <- list.files(dir,pattern="HepG2",ignore.case = T)
  combine_bw(dir,flist_hepg2,paste0(out_dir,out_prefix,"_HepG2.bw"),"avg",num_samples)
  flist_k562 <- list.files(dir,pattern="K562",ignore.case = T)
  combine_bw(dir,flist_k562,paste0(out_dir,out_prefix,"_K562.bw"),"avg",num_samples)
}

process_bw_wrapper("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/DHS/BW","DHS")
process_bw_wrapper("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/BW_hg38/","SuRE")


