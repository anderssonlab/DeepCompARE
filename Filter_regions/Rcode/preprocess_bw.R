source("/binf-isilon/alab/people/pcr980/R/Function_bigwig.R")

#------------------------------------
# average bigwig files of replicates
#-----------------------------------
out_dir <- "/isdata/alab/people/pcr980/Raw_data/BW_processed/"

process_bw_wrapper <- function(dir,out_prefix,num_samples="infer"){
  flist_hepg2 <- list.files(dir,pattern="HepG2",ignore.case = T)
  combine_bw(dir,flist_hepg2,paste0(out_dir,out_prefix,"_HepG2.bw"),"avg",num_samples)
  flist_k562 <- list.files(dir,pattern="K562",ignore.case = T)
  combine_bw(dir,flist_k562,paste0(out_dir,out_prefix,"_K562.bw"),"avg",num_samples)
}

process_bw_wrapper("/isdata/alab/people/pcr980/Raw_data/DHS_files/BW_files/","DHS")
process_bw_wrapper("/isdata/alab/people/pcr980/Raw_data/CAGE_files/BW_files_TPM/","CAGE",num_samples=6)
process_bw_wrapper("/isdata/alab/people/pcr980/Raw_data/SuRE_files/BW_files/","SuRE")


