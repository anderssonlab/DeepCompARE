setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions")
library(rtracklayer)
source("Rcode/Function_bigwig.R")
    # combine_bw

process_dir <- function(dir,cell_type,outname){
        files <- list.files(dir)
        files <- grep(cell_type,files,value = TRUE)
        combine_bw(dir,files,outname,"sum")
}

process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE42_HG02601",
            "HEPG2",
            "SuRE/BW_hg19/sure42_hepg2_hg19.bw")
process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE42_HG02601",
            "K562",
            "SuRE/BW_hg19/sure42_k562_hg19.bw")

process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE43_GM18983",
            "HEPG2",
            "SuRE/BW_hg19/sure43_hepg2_hg19.bw")
process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE43_GM18983",
            "K562",
            "SuRE/BW_hg19/sure43_k562_hg19.bw")

process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE44_HG01241",
            "HEPG2",
            "SuRE/BW_hg19/sure44_hepg2_hg19.bw")
process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE44_HG01241",
            "K562",
            "SuRE/BW_hg19/sure44_k562_hg19.bw")

process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE45_HG03464",
            "HEPG2",
            "SuRE/BW_hg19/sure45_hepg2_hg19.bw")
process_dir("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/Signal_SuRE45_HG03464",
            "K562",
            "SuRE/BW_hg19/sure45_k562_hg19.bw")


