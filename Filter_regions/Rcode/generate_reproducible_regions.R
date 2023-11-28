setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions")
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(patchwork)

source("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/Rcode/Function_EDA.R")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/Rcode/Function_utils.R")
source("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/Rcode/Function_reproducibility.R")

modalities <- c("CAGE/","DHS/","STARR/","SuRE/") # The order is always like this.
cells <- c("HepG2","K562") # The order is always like this.

dir_base <- "/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/"
dir_EDA <- paste0(dir_base,"EDAs/")

dirs<- c(paste0(dir_base,modalities))
dirs_reproducible <- paste0(dirs,"Peaks_reproducible_within_modality")

dirs_irreproducible <- paste0(dirs,c("Peaks_irreproducible",
                                     "Peaks_irreproducible",
                                     "Peaks_irreproducible",
                                     "NarrowPeak_hg19_c4_RMSNP"))


# No need to generate loose for CAGE and STARR, already generated in their separate files.
# Generate 
dirs_loose <- paste0(dirs,c("","Peaks_irreproducible","","NarrowPeak_hg38_c2"))


#------------------------------------------
# Step 0: EDA on original file
#------------------------------------------
# check validity of original file
# output region width, dist2nearest within file, dist2nearest rDHS(ENCFF619EJE.bed), intersection/union rDHS

dirs_irreproducible[4] <- paste0(dirs[4],"NarrowPeak_hg38_c4_RMSNP") # use hg38 to compare with rDHS
files <- list.files(dirs_irreproducible,full.name=T)
EDA_path <- paste0(dir_EDA,"1.validity_single_original_file.pdf")
pdf(EDA_path,height = 7, width = 28)
lapply(files,single_file_validity_EDA)
dev.off()



# reproducibility within modality of original files
dirs_irreproducible[4] <- paste0(dirs[4],"NarrowPeak_hg19_c4_RMSNP") # use hg19 to plot intersection lengths
EDA_path <- paste0(dir_EDA,"2.reproducibility_within_modality_of_original_files.pdf")
pdf(EDA_path,height = 4, width = 8)
reproducibility_within_modality_EDA(dirs_irreproducible[2],120,"DHS")
reproducibility_within_modality_EDA(dirs_irreproducible[3],400,"STARR")
reproducibility_within_modality_EDA(dirs_irreproducible[4],240,"SuRE")
dev.off()



#-----------------------------------------------------
# Step1: output regions reproducible within modality
# SuRE data should use hg19 version to avoid the effect of liftOver
#-----------------------------------------------------
dirs_irreproducible[4] <- paste0(dirs[4],"NarrowPeak_hg19_c4_RMSNP") # use hg19 to output intersection as reproducible regions

# DHS
export_reproducible_within_modality(dirs_irreproducible[2],120, 2, 
                                    paste0(dirs_reproducible[2],"/DHS"))
# STARR
export_reproducible_within_modality(dirs_irreproducible[3],400, 2, 
                                    paste0(dirs_reproducible[3],"/STARR"))
# SuRE: First run macs3 bdgpeakcall -c 4, then SuRE_narrowPeak_EDA_filter_liftOver.R
export_reproducible_within_modality(dirs_irreproducible[4],240,2,
                                     paste0(dirs_reproducible[4],"/SuRE"),lift=T)



EDA_path <- paste0(dir_EDA,"3.validity_single_modality_reproducible_regions_within_modality.pdf")
pdf(EDA_path,height = 7, width = 28)
lapply(dirs_reproducible,single_modality_validity_EDA)
dev.off()

EDA_path <- paste0(dir_EDA,"4.modality_compatibility_for_regions_within_modality.pdf")
pdf(EDA_path,height = 7, width = 28)
modality_compatibility_EDA(dirs_reproducible)
dev.off()


#-------------------------------------
# Step3: output rescued peaks
# regions reproducible within modality format: bed
# regions irreproducible format: narrowPeak
# SuRE data should use hg38 version
#-------------------------------------
dirs_irreproducible[4] <- paste0(dirs[4],"NarrowPeak_hg38_c4_RMSNP") # use hg38 to rescue other modalities
modalities <- c("CAGE","DHS","STARR","SuRE")
for (i in 1:4){
  gr_hepg2 <- single_celltype_supported_across_modality(dirs_irreproducible[i],dirs_reproducible[i],"HepG2")
  export.bed(gr_hepg2,paste0(dirs[i],"Peaks_rescued/",modalities[i],"_HepG2.bed"))
  
  gr_k562 <- single_celltype_supported_across_modality(dirs_irreproducible[i],dirs_reproducible[i],"K562")
  export.bed(gr_k562,paste0(dirs[i],"Peaks_rescued/",modalities[i],"_K562.bed"))
}


EDA_path <- paste0(dir_EDA,"5.validity_single_file_rescued_regions.pdf")
pdf(EDA_path,height = 7, width = 28)
lapply(list.files(paste0(dirs,"Peaks_rescued/"),full.name=T),single_file_validity_EDA)
dev.off()


EDA_path <- paste0(dir_EDA,"6.validity_single_modality_rescued_regions.pdf")
pdf(EDA_path,height = 7, width = 28)
lapply(paste0(dirs,"Peaks_rescued/"),single_modality_validity_EDA)
dev.off()

EDA_path <- paste0(dir_EDA,"7.modality_compatibility_for_rescued_regions.pdf")
pdf(EDA_path,height = 7, width = 28)
modality_compatibility_EDA(paste0(dirs,"Peaks_rescued/"))
dev.off()


#---------------------------------------------------------------
# Step3: combine peaks reproducible within modality and rescued
# Res dir: /maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/
#---------------------------------------------------------------


for (i in 1:4){
  dir <- dirs[i]
  lapply(cells,function(cell){
    out_name <- paste0("/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/",modalities[i],"_",cell,".bed")
    combine_peaks_reproducible_and_rescued(dir,cell,out_name)
  })
}


#----------------------------------------------------------------------
# Step 4: combine loose positive for each modality
#----------------------------------------------------------------------
export_loose_within_modality(dirs_loose[2],"DHS")
export_loose_within_modality(dirs_loose[4],"SuRE")
