library(dplyr)

# plot feature attribution
plot_attr <- function(position,feat_attr){
  feat_attr <- as.numeric(feat_attr)
  position <- as.numeric(position)
  df_plot <- as.data.frame(cbind(position,feat_attr))
  colnames(df_plot) <- c("position","attribution")
  ggplot(df_plot,aes(x=position, y=attribution)) +
    geom_bar(fill="darkgrey", stat="identity",position=position_identity(), width=.8) + 
    xlab("") + theme_bw()
}

# experimental value
df <- read.csv("Validation_SNP_MPRA/MPRA_Enhancer_HepG2_SORT1.csv")
my_group <- paste0(df$Ref,df$relevant_position)
df <- aggregate(df$Value, by = list(my_group), FUN = mean)
df$position <- as.numeric(gsub("\\D", "", df$Group.1))
ggplot(df,aes(position, x)) +
  geom_bar(fill="darkgrey", stat="identity",position=position_identity(), width=.8) + 
  xlab("") + theme_bw()








#-----------------
# Plot gradientxInput of SORT1 enhancer
#----------------
df_grad <- read.csv("Showcase_SORT1/gradxinp_track_starr_hepg2.csv",
                    row.names = 1,
                    header=F)
rownames(df_grad) <- NULL
df_seq <- read.csv("Showcase_SORT1/seq_info.csv")
sequence_ref <- unlist(strsplit(df_seq$seq[1], split = ""))
sequence_alt <- unlist(strsplit(df_seq$seq[2], split = ""))

df_plot1 <- as.data.frame(cbind(200:400,as.numeric(df_grad[1,200:400]),sequence_ref[200:400],"SORT1_enhancer"))
colnames(df_plot1) <- c("Position","Importance","Base","Sequence")
df_plot2 <- as.data.frame(cbind(200:400,as.numeric(df_grad[2,200:400]),sequence_alt[200:400],"Mutated_SORT1_enhancer"))
colnames(df_plot2) <- c("Position","Importance","Base","Sequence")

df_plot <- as.data.frame(rbind(df_plot1,df_plot2))
df_plot$Position <- as.numeric(df_plot$Position)
df_plot$Importance <- as.numeric(df_plot$Importance)

mutation_data <- data.frame(
  Position = 317,
  Importance = as.numeric(df_plot2$Importance[df_plot2$Position == 317]),
  Sequence = "Mutated_SORT1_enhancer"
)


ggplot(df_plot,aes(x=Position, y=Importance)) +
  geom_bar(fill="darkgrey", stat="identity",position=position_identity(), width=.8) + 
  geom_point(aes(color=Base),size=0.5) +
  facet_wrap(~Sequence,ncol=1,scales="free_y") +
  geom_rect(aes(xmin = 316, xmax = 324, ymin = -Inf, ymax = Inf),fill = "pink",alpha = 0.01)+
  geom_point(data = mutation_data, aes(x = Position, y = Importance),
             color = "red", size = 2) +
  geom_text(data = mutation_data, aes(x = Position, y = Importance),
            label = "rs12740374:G>T", vjust = -1.2,hjust=1, color = "red", size = 4) +
  xlab("Relative position") + 
  ylab("GradxInp from STARR HepG2 track")+
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size=4)))










