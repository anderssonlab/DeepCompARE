setwd("/maps/projects/ralab/people/pcr980/DeepCompare/SNP_interpretation_prime_editing")
library(ggplot2)


df <- read.csv("correlation.csv")
colnames(df) <- c("Method","Correlation","Pval", "Regulatory_element", "Model")

df <- df[!grepl("hepg2",df$Method),]

df$Method <- gsub("_regression","",df$Method)
df$Method <- gsub("_classification","",df$Method)
df$Method <- gsub("_k562","",df$Method)

df$Method <- gsub("CAGE","cage",df$Method)
df$Method <- gsub("DNase","dhs",df$Method)

df$Model <- gsub("DeepCompare regression","DeepCompare\nregression",df$Model)
df$Model <- gsub("DeepCompare classification","DeepCompare\nclassification",df$Model)


p <- ggplot(df,(aes(x=Model,y=Correlation,color=Model,shape=Method)))+
  geom_point(position = position_dodge(width = 0.5),size=3)+
  scale_shape_manual(values = c(16,17,4,8))+
  facet_wrap(~ Regulatory_element, scales = "free_y", ncol = 2)+
  theme(axis.text.x = element_blank(), 
        axis.text.x.bottom = element_text(),
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove panel background
        plot.background = element_blank(),
        #panel.grid.major = element_line(color = "grey"),
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 10)
  ) +
  labs(title = "Tile prime editing: DeepCompare v.s. Enformer", 
       x = "Model", 
       y = "Pearson Correlation")


ggsave(paste0("Tile_prime_editing_DeepCompare_vs_Enformer.pdf"), plot = p, width = 10, height = 4)
