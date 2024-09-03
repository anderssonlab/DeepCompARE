setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed")
library(rtracklayer)
library(GenomicRanges)

files_exclude <- list.files(pattern="_loose.bed")
files_exclude <- c(files_exclude,
                   "hg38_non_hotspots_non_blacklist_600bp_windows_intergenic.bed",
                   "hg38_non_hotspots_non_blacklist_600bp_windows.bed",
                   "resize.R",
                   "universe.bed",
                   "SuRE_SNP_hg38.bed")


files_include <- setdiff(list.files(),files_exclude)

lapply(files_include,function(file){
  bed <- import.bed(file)
  bed <- resize(bed,width=600,fix = "center")
  export(bed, paste0("resize_600bp_",file))
})
