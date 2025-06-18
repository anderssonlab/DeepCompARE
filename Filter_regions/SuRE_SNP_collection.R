setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Filter_regions/")
library(rtracklayer)

#-------------------------------------
# Summarize SNP data for each sample
#-------------------------------------

get_snp <- function(dir_name,fname){
	df <- read.table(paste0(dir_name,fname),sep="\t",header=TRUE)
	df <- df[!df$SNPrelpos=="",]
	abs_pos <- unique(df$SNPabspos)
	abs_pos <- as.numeric(unlist(strsplit(abs_pos,split=',', fixed=TRUE)))
	gr <- GRanges(seqnames = unique(df$chr),ranges = IRanges(abs_pos))
	unique(gr)
}

report_snps <- function(dir_name,output_name){
	files <- list.files(dir_name)
	res <- sapply(files, function(fname){get_snp(dir_name,fname)})
	res_flattened <- res[[1]]
	for (i in 2:length(res)){
		res_flattened <- c(res_flattened,res[[i]])
	}
	res_flattened <- unique(res_flattened)
	export.bed(res_flattened,con=output_name)
}

report_snps("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/SNP_SuRE42_HG02601/","SuRE/SuRE_SNP/SNP_SuRE42_hg19.bed")
report_snps("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/SNP_SuRE43_GM18983/","SuRE/SuRE_SNP/SNP_SuRE43_hg19.bed")
report_snps("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/SNP_SuRE44_HG01241/","SuRE/SuRE_SNP/SNP_SuRE44_hg19.bed")
report_snps("/maps/projects/ralab/people/pcr980/DeepCompare/Raw_data/SuRE/SNP_SuRE45_HG03464/","SuRE/SuRE_SNP/SNP_SuRE45_hg19.bed")


#-----------------------------------
# Combine SNP info of all samples
#-----------------------------------


dir <- "SuRE/SuRE_SNP/"

f42 <- import(paste0(dir,"SNP_SuRE42_hg38.bed"),format="BED")
f43 <- import(paste0(dir,"SNP_SuRE43_hg38.bed"),format="BED")
f44 <- import(paste0(dir,"SNP_SuRE44_hg38.bed"),format="BED")
f45 <- import(paste0(dir,"SNP_SuRE45_hg38.bed"),format="BED")
gr <- c(f42,f43,f44,f45)
gr <- unique(gr)
mcols(gr) <- NULL
export.bed(gr,"/maps/projects/ralab/people/pcr980/DeepCompare/Pd1_bed_processed/SuRE_SNP_hg38.bed")



