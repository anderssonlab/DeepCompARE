setwd("/maps/projects/ralab/people/pcr980/DeepCompare/RNA_expression/")

library(rtracklayer)
library(biomaRt)
library(dplyr)

# read in expression data
expr_hepg2 <- read.table("Raw_data/ENCFF103FSL_RNASeq_HepG2.tsv",header = T)
expr_k562 <- read.table("Raw_data/ENCFF928NYA_RNASeq_K562.tsv",header = T)

expr_hepg2$ensembl <- ifelse(grepl("^ENSG", expr_hepg2$gene_id), sub("\\.[0-9]+$", "", expr_hepg2$gene_id), NA)
expr_k562$ensembl <- ifelse(grepl("^ENSG", expr_k562$gene_id), sub("\\.[0-9]+$", "", expr_k562$gene_id), NA)


# read in TF-gene name mapper
name_mapper <- read.csv("/maps/projects/ralab/people/pcr980/Resource/HTF/DatabaseExtract_v_1.01.csv",row.names = 1)
sum(name_mapper$Is.TF.=="Yes")
name_mapper <- name_mapper %>% filter(grepl("^ENSG",Ensembl.ID))
rownames(name_mapper) <- name_mapper$Ensembl.ID



# subset expression data by HTF
expr_hepg2 <- expr_hepg2 %>% filter(ensembl %in% name_mapper$Ensembl.ID)
expr_k562 <- expr_k562 %>% filter(ensembl %in% name_mapper$Ensembl.ID)
rownames(expr_hepg2) <- expr_hepg2$ensembl
rownames(expr_k562) <- expr_k562$ensembl
all(rownames(expr_hepg2)==rownames(expr_k562))
expr_hepg2$ensembl <- NULL
expr_k562$ensembl <- NULL

name_mapper <- name_mapper[rownames(expr_hepg2),]
all(rownames(name_mapper)==rownames(expr_hepg2))

tf_expr_hepg2 <- cbind(expr_hepg2,name_mapper)
tf_expr_k562 <- cbind(expr_k562,name_mapper)

all(rownames(tf_expr_hepg2)==tf_expr_hepg2$Ensembl.ID)
all(rownames(tf_expr_k562)==tf_expr_k562$Ensembl.ID)

write.csv(tf_expr_hepg2,"TF_expression_hepg2.csv")
write.csv(tf_expr_k562,"TF_expression_k562.csv")




#-------------------------------------------------------
# Validation: Are TF with ENCODE ChIP experiments expressed?
#-------------------------------------------------------

# remap: TF with ENCODE ChIP experiments in HepG2 and K562
remap_hepg2 <- import.bed("/maps/projects/ralab/people/pcr980/Resource/ReMap2022/ReMap2022_Hep-G2.bed")
remap_k562 <- import.bed("/maps/projects/ralab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")

tf_hepg2 <-  sub(":.*", "", remap_hepg2$name)
tf_hepg2 <- unique(tf_hepg2)
tf_k562 <-  sub(":.*", "", remap_k562$name)
tf_k562 <- unique(tf_k562)

tf_expr_hepg2_sub <- tf_expr_hepg2 %>% filter(HGNC.symbol %in% tf_hepg2) # two dropout
tf_expr_k562_sub <- tf_expr_k562 %>% filter(HGNC.symbol %in% tf_k562) # two dropout



#-------------------------------------------------------
# Create a list of expressed and nonexpressed TFs for hepG2 and K562
#-------------------------------------------------------
tf_expressed_hepg2 <- tf_expr_hepg2 %>% filter(TPM>1) %>% pull(HGNC.symbol) 
tf_expressed_hepg2 <- unique(c(tf_expressed_hepg2,tf_hepg2))

tf_expressed_k562 <- tf_expr_k562 %>% filter(TPM>1) %>% pull(HGNC.symbol) 
tf_expressed_k562 <- unique(c(tf_expressed_k562,tf_k562))


tf_unexpressed_hepg2 <- tf_expr_hepg2 %>% filter(TPM==0) %>% pull(HGNC.symbol) 
tf_unexpressed_hepg2 <- tf_unexpressed_hepg2[!tf_unexpressed_hepg2 %in% tf_hepg2]

tf_unexpressed_k562 <- tf_expr_k562 %>% filter(TPM==0) %>% pull(HGNC.symbol) 
tf_unexpressed_k562 <- tf_unexpressed_k562[!tf_unexpressed_k562 %in% tf_k562]

write.table(tf_expressed_hepg2,"expressed_tf_list_hepg2.tsv",row.names = F,col.names = F)
write.table(tf_expressed_k562,"expressed_tf_list_k562.tsv",row.names = F,col.names = F)

write.table(tf_unexpressed_hepg2,"unexpressed_tf_list_hepg2.tsv",row.names = F,col.names = F)
write.table(tf_unexpressed_k562,"unexpressed_tf_list_k562.tsv",row.names = F,col.names = F)
