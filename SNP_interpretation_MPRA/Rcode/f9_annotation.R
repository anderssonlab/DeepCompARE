setwd("~/Documents/DeepCompare/SNP_interpretation_MPRA")
library(Biostrings)


# USCS genome annotation downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
gtf <- read.table('hg38.knownGene.gtf', header = FALSE, sep = '\t')

# F9 position: chrX:139,530,463-139,530,765, original width 303
# F9 gene_name P00740
# F9 TSS: 139,530,739, from GTF file
# F9 +1 relevant position: 423
# F9 -1 relevant position: 422

df <- read.csv("Pd1_MPRA_with_seq_info/ref_seq_forward.csv",row.names = 1)
ref_seq <- df[df$RE_name=="F9","ref_seq"]

ref_seq <- DNAString(ref_seq)


# CEBPA: +6, relevant position: 423+6=429
ref_seq[428:435]

# ONECUT: -9, relevant position: 422-8=414
ref_seq[414:419]

# HNF4A: -27, relevant position: 422-26=396
ref_seq[396:404]

# AR: -36, relevant position: 422-35=387
ref_seq[387:401]
