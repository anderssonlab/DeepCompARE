library(ggplot2)
# Read data
# Only investigate CAGE HepG2
df_classification <-read.csv("MPRA_paper_deepCompare_annotated_results_CLASSIFICATION.csv")
df_regression <- read.csv("MPRA_paper_deepCompare_annotated_results_REGRESSION.csv")
df_classification <- df_classification[df_classification$celltype=="hepg2",]
df_classification <- df_classification[df_classification$assay=="cage",]
df_regression <- df_regression[df_regression$celltype=="hepg2",]
df_regression <- df_regression[df_regression$assay=="cage",]

# filter out low tags
df_classification <- df_classification[df_classification$tags>10,]
df_regression <- df_regression[df_regression$tags>10,]

# filter out scrambled
idx <- which(df_classification$Scrambled==TRUE)
df_classification <- df_classification[-idx,]
df_regression <- df_regression[-idx,]

# filter out "other"
idx <- which(df_classification$type=="other")
df_classification <- df_classification[-idx,]
df_regression <- df_regression[-idx,]

# filter out construction 2 (lower signal to noise ratio):
idx <- which(df_classification$construct.x=="Construct2")
df_classification <- df_classification[-idx,]
df_regression <- df_regression[-idx,]


cor(df_classification$score_dc,df_classification$ratio)
cor(df_regression$score_dc,df_regression$ratio)

df_sub <- df_classification[df_classification$type=="Single Motif",]
cor(df_sub$score_dc,df_sub$ratio)

# plot
ggplot(data=df_classification,aes(x=log2(ratio),y=score_dc,color=as.factor(HNF1A)))+
  geom_point(size=0.1,alpha=0.1)+
  guides(color=guide_legend(override.aes=list(size=4, alpha=1)))+
  theme_bw()
  


df_subsub<- df_sub[df_sub$type=="Single Motif",]
cor(df_subsub$score_dc,df_subsub$ratio)

# See the nike
cor(df_sub$score_dc,df_sub$ratio)
plot(df_classification$score_dc,df_regression$score_dc)



# Homo is enriched in short arm (high regression low classification)
homo_idx <- grep("Homo",df_classification$type) # 20314 homo/191554 total
idx <- which(df_classification$score_dc<10 & df_regression$score_dc>1) # select the short arm
df_short <- df_classification[idx,]
homo_short_idx <- grep("Homo",df_short$type) # 477 homo/1193 total



# Homos in long arm don't have REST, YY1, GABPA, CTCF, ONECUT, AHR, TFAP2C, CREB1,AP1, SP1,XBP1, NR2F2,PPARA
df_homo_class <- df_classification[homo_idx,]
df_homo_reg <- df_regression[homo_idx,]
plot(df_homo_class[,"score_dc"],df_homo_reg[,"score_dc"]) # long arm gets thinner


idx <- which(df_homo_class$score_dc>10) # select the long arm
df_homo_long <- df_homo_class[idx,]
# For homo, except the following listed, all the rest TFBS are absent in the long arm
table(df_homo_long$CEBPA) # CEBPA is almost absent
table(df_homo_long$RXRA) 
table(df_homo_long$HNF4A) 
table(df_homo_long$FOXA1) 
table(df_homo_long$HNF1A) 



# Hetero in long arm have maximum 3 REST, YY1, GABPA, CTCF, ONECUT, AHR, TFAP2C, CREB1, AP1, XPBP1, CEBPA, NR2F2, PRARA (Inhibitory?)
# Hetero in long arm have maximum 2 SP1
# Hetero in long arm have maximum 6 RXRA, HNF4A, FOXA1, HNF1A (Pioneering excitatory?)
df_class_hetero <- df_classification[-homo_idx,]
idx <- which(df_class_hetero$score_dc>10)
df_hetero_long <- df_class_hetero[idx,]
table(df_hetero_long$Scrambled)

df_hetero_long$excitatory <- df_hetero_long$HNF1A+df_hetero_long$HNF4A+df_hetero_long$RXRA+df_hetero_long$FOXA1

cor(df_hetero_long$score_dc,df_hetero_long$HNF1A)
cor(df_hetero_long$score_dc,df_hetero_long$HNF4A)
cor(df_hetero_long$score_dc,df_hetero_long$RXRA)
cor(df_hetero_long$score_dc,df_hetero_long$FOXA1)
cor(df_hetero_long$score_dc,df_hetero_long$excitatory)

plot(df_hetero_long$REST,df_hetero_long$excitatory)

plot(df_hetero_long$score_dc,df_hetero_long$ratio)
