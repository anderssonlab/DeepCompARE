setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Test")
library(ggplot2)
library(dplyr)




#------------------
# ST models
#-----------------

df_st <- read.csv("ST_metrics.csv")
df_st$X <- NULL

# Compare Class and CR for Conv5D
df_sub <- df %>% 
  filter(model_name %in% c("Conv5D","CRConv5D")) %>% 
  filter(model_type %in% c("Class", "CR"))

ggplot(df_sub,aes(x=file,y=acc,fill=model_type))+
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "ACC", x = "File") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

# Compare Reg and CR for Conv5D
df_sub <- df %>% 
  filter(model_name %in% c("Conv5D","CRConv5D")) %>% 
  filter(model_type %in% c("Reg", "CR"))

ggplot(df_sub,aes(x=file,y=pcc,fill=model_type))+
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "PCC", x = "File") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")





#------------------
# MT models
#-----------------
df_mt <- read.csv("MT_metrics.csv",header=FALSE)
colnames(df_mt) <- c('file', 'pcc', 'acc', 'model_name', 'task_type', 'model_type', 'data_dir')

df_mt$data_dir <- gsub("_[0-9]$","",df_mt$data_dir)
df_mt <- df_mt %>% group_by(data_dir,model_type,file) %>% 
  summarize(mean_pcc=mean(pcc),mean_acc=mean(acc),sd_pcc=sd(pcc),sd_acc=sd(acc))

# compare with and without multilabel
df_mt_sub <- df_mt %>% filter(data_dir!="Dataset_5cv_loose" & model_type=="CR")
df_mt_sub$data_dir <- gsub("Dataset_5cv_without_multilabel","Dataset without multilabel",df_mt_sub$data_dir)
df_mt_sub$data_dir <- gsub("Dataset_5cv_with_multilabel","Dataset with multilabel",df_mt_sub$data_dir)
colnames(df_mt_sub)[1] <- "Training_dataset"
p <- ggplot(df_mt_sub,aes(x=file,y=mean_pcc,color=Training_dataset))+
  geom_point(position = position_dodge(width = 0.5),size=1)+
  geom_errorbar(aes(ymin = mean_pcc-sd_pcc, ymax = mean_pcc+sd_pcc),position = position_dodge(width = 0.5),width = 0)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Dataset")+
  ylab("Pearson correlation")
ggsave("Generated_plots/With_multilabel_is_better.pdf", plot = p, width = 8, height = 4)



# remove dataset without multilabel
df_mt <- df_mt[df_mt$data_dir!="Dataset_5cv_without_multilabel",]
colnames(df_mt)[2] <- "Model" 
df_mt$Model <- gsub("Reg","Regression",df_mt$Model)
df_mt$Model <- gsub("Class","Classification",df_mt$Model)
df_mt$Model <- factor(df_mt$Model,levels=c("Regression","Classification","CR"))

p <- ggplot(df_mt[df_mt$data_dir=="Dataset_5cv_with_multilabel" & df_mt$Model %in% c("CR","Regression"),], 
       aes(x=file,y=mean_pcc,color=Model))+
  geom_point(position = position_dodge(width = 0.5),size=1)+
  geom_errorbar(aes(ymin = mean_pcc-sd_pcc, ymax = mean_pcc+sd_pcc),position = position_dodge(width = 0.5),width = 0)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Dataset")+
  ylab("Pearson correlation")
ggsave("Generated_plots/CR_better_than_reg.pdf", plot = p, width = 8, height = 4)


p <- ggplot(df_mt[df_mt$data_dir=="Dataset_5cv_with_multilabel" & df_mt$Model %in% c("CR","Classification"),], 
       aes(x=file,y=mean_acc,color=Model))+
  geom_point(position = position_dodge(width = 0.5),size=1)+
  geom_errorbar(aes(ymin = mean_acc-sd_acc, ymax = mean_acc+sd_acc),position = position_dodge(width = 0.5),width = 0)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Dataset")+
  ylab("Accuracy")
ggsave("Generated_plots/CR_comparable_to_classification.pdf", plot = p, width = 8, height = 4)


# Output: compare CR-reg branch vs loose
df_mt <- df_mt[df_mt$data_dir=="Dataset_5cv_loose" | df_mt$Model=="CR",]
df_mt$Model <- as.character(df_mt$Model)
df_mt$Model <- gsub("Regression","Regression trained on all regions",df_mt$Model)
df_mt$Model <- gsub("CR","CR trained on credible regions",df_mt$Model)
df_mt$Model <- factor(df_mt$Model,levels=c("Regression trained on all regions","CR trained on credible regions"))
p <- ggplot(df_mt,aes(x=file,y=mean_pcc,color=Model))+
  geom_point(position = position_dodge(width = 0.5),size=1)+
  geom_errorbar(aes(ymin = mean_pcc-sd_pcc, ymax = mean_pcc+sd_pcc),position = position_dodge(width = 0.5),width = 0)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Dataset")+
  ylab("Pearson correlation")

  ggsave("Generated_plots/Credible_regions_outperform_loose.pdf", plot = p, width = 8, height = 4)


