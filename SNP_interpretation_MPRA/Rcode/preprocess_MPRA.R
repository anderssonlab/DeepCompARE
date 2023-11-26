#Add reference and alternative sequence (both forward and reverse strand) to data set

#Elements of interest has following location:

#F9:chrX:139,530,463-139,530,765
#LDLR:chr19:11,089,231-11,089,548
#PKLR: chr1:155,301,395-155,301,864
#SORT1: chr1:109,274,652-109,275,251

# GP1BA chr22:19,723,266-19,723,650
# HBB chr11:5,227,022-5,227,208
# HBG1 chr11:5,249,805-5,250,078
# BCL11A	chr2:60,494,940-60,495,539

setwd("~/Documents/DeepCompare/SNP_interpretation_MPRA") # if working on local computer
#setwd("/binf-isilon/alab/people/pcr980/SNP_interpretation_MPRA") # if working on binf server
source("~/Documents/DeepCompare/Rcode/variant_interpretation.R")
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

# read in raw data
df <- read.csv("Raw_data/GRCh38_F9_GP1BA_HBB_HBG1_LDLR_LDLR.2_PKLR-24h_PKLR-48h_BCL11A_SORT1_SORT1-flip_SORT1.2.csv")
dim(df) # 17489 variant before preprocessing
unique(df$Element)
df$RE_name <- gsub("\\..*|\\-.*", "", df$Element)


# get reference sequence (forward strand)
df_gr <- data.frame(
  RE_name=c("F9","LDLR","PKLR","SORT1","GP1BA","HBB","HBG1","BCL11A"),
  seqnames=c("chrX","chr19","chr1","chr1","chr22","chr11","chr11","chr2"),
  start=c(139530463,11089231,155301395,109274652,19723266,5227022,5249805,60494940),
  end=c(139530765,11089548,155301864,109275251,19723650,5227208,5250078,60495539)
)
df_gr$length <- df_gr$end-df_gr$start+1
gr <- makeGRangesFromDataFrame(df_gr)
gr <- resize(gr,600,fix="center")
strand(gr) <- "+"
df_gr$ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr)


# write down reference strand info, for Enformer
write.csv(df_gr,"Pd1_MPRA_with_seq_info/ref_seq_forward.csv")
strand(gr) <- "-"
df_gr$ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr)
write.csv(df_gr,"Pd1_MPRA_with_seq_info/ref_seq_reverse.csv")

# merge 
df <- merge(df,df_gr,by="RE_name")
df$seqnames <- NULL

# filter out low-quality rows. After filtering 5860 rows remain.
df <- df[df$P.Value<0.01,]
df <- df[df$Tags>10,]
df <- df[df$Alt!="-",]

# get relevant position of mutation (forward strand)
df$relevant_position <- df$Position-df$start+1+ceiling((600-df$length)/2)
all(sapply(1:dim(df)[1],function(i) df[i,"Ref"]==as.character(df[[i,"ref_seq"]][df[i,"relevant_position"]])))

df$alt_seq <- sapply(1:dim(df)[1],function(i){
  this_ref_seq <- df[i,"ref_seq"]
  relevant_pos <- df[i,"relevant_position"]
  this_alt_base <- df[i,"Alt"]
  as.character(replaceAt(this_ref_seq,IRanges(relevant_pos,relevant_pos),this_alt_base))
})


# write down forward strand info
write.csv(df,"Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv")


# write down reverse strand info

df_rc <- df
df_rc$ref_seq <- get_rc_from_char_vec(df_rc$ref_seq)
df_rc$alt_seq <- get_rc_from_char_vec(df_rc$alt_seq)
write.csv(df_rc,"Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv")




