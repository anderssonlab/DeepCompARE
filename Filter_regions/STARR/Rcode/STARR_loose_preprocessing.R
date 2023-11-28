setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/Rcode/Function_utils.R")


data_dir <- "STARR/STARRPeaker_res/"

gr <- read_file(paste0(data_dir,"STARR_HepG2.peak.final.starrPeak"))
gr <- gr[gr$pval>-log10(0.2)] 
export.bed(gr,"/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/STARR_HepG2_loose.bed")
rm(gr)

gr <- read_file(paste0(data_dir,"STARR_K562.peak.final.starrPeak"))
gr <- gr[gr$pval>-log10(0.2)] 
export.bed(gr,"/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/STARR_K562_loose.bed")
