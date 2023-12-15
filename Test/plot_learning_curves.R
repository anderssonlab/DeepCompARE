setwd("/maps/projects/ralab/people/pcr980/DeepCompare/Test")
library(dplyr)
library(ggplot2)
#--------------------
# Functions
#--------------------
split <- function(model_name){
  words <- unlist(strsplit(model_name,"_"))
  model_type <- words[length(words)-1]
  dataset <- paste(words[2:(length(words)-3)],collapse="_")
  return(c(model_type,dataset))
}


process_one_model <- function(model){
  df_train<- read.csv(file.path(dir_models,model,"lc_train_loss.csv"),row.names = 1)
  df_train$loss <- rowMeans(df_train)
  df_train$epoch <- 1:100
  df_train$type <- "train"
  df_train_acc <- read.csv(file.path(dir_models,model,"lc_train_acc.csv"),row.names = 1)
  df_train$acc <- df_train_acc[,1]
  
  df_val <- read.csv(file.path(dir_models,model,"lc_val_loss.csv"),row.names = 1)
  df_val$loss <- rowMeans(df_val)
  df_val$epoch <- seq(2,100,2)
  df_val$type <- "val"
  df_val_acc <- read.csv(file.path(dir_models,model,"lc_val_acc.csv"),row.names = 1)
  df_val$acc <- df_val_acc[,1]
  
  df <- rbind(df_train,df_val)
  df <- df[,c("loss","acc","epoch","type")]
  info <- split(model)
  df$model_type <- info[1]
  df$dataset <- info[2]
  
  df
}

#---------------------------------------
# Aggregate information for all models
#---------------------------------------

dir_models <- "/maps/projects/ralab/people/pcr980/DeepCompare/Models"
models <- list.files(dir_models)
models <- models[grepl("Astig",models)]

dfs <- lapply(models,process_one_model)
df <- do.call(rbind,dfs)


#---------------------
# Plot
#---------------------
df$dataset <- gsub("_[0-9]$","",df$dataset)
df_plot <- df %>% 
  group_by(dataset,type,model_type,epoch) %>% 
  summarize(median_loss=median(loss),sd_loss=sd(loss),
            median_acc=median(acc),sd_acc=sd(acc))

# with/without multilabel
df_sub <- df_plot%>%filter(dataset!="Dataset_5cv_loose" & model_type=="CR")
ggplot(df_sub, aes(x=epoch,y=median_acc,color=dataset,linetype=type))+
  geom_line()

# Reg v.s. CR
df_sub <- df_plot%>%filter(dataset=="Dataset_5cv_with_multilabel" & model_type %in% c("CR","Reg"))
ggplot(df_sub, aes(x=epoch,y=median_loss,color=model_type,linetype=type))+
  geom_line()


# Class v.s. CR
df_sub <- df_plot%>%filter(dataset=="Dataset_5cv_with_multilabel" & model_type %in% c("CR","Class"))
ggplot(df_sub, aes(x=epoch,y=median_acc,color=model_type,linetype=type))+
  geom_line()
