library(ggplot2)

plot_importance_by_position_single_track <- function(positions,importances,bases=NULL){
  # plot importance by position for single track of importance
  positions <- as.numeric(positions)
  importances <- as.numeric(importances)
  if (is.null(bases)){
    df_plot <- data.frame(list(Position=positions,Importance=importances))
    ggplot(df_plot,aes(x=Position, y=Importance)) +
      geom_bar(fill="darkgrey", stat="identity",position=position_identity(), width=1) + 
      xlab("Relevant position") + 
      ylab("Effect size")+
      geom_hline(yintercept = 0, color = "black", linewidth = 1) +
      theme_minimal()+
      theme(
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()
      )
  } else{
    bases <- as.character(bases)
    df_plot <- data.frame(list(Position=positions,Importance=importances,Alternative_base=bases))
    ggplot(df_plot,aes(x=Position, y=Importance)) +
      geom_bar(fill="darkgrey", stat="identity",position=position_identity(), width=1) + 
      geom_point(df_plot,aes(x=relevant_position, y=Importance, color=Alternative_base),size=0.6)+
      guides(color = guide_legend(override.aes = list(size=2)))+
      xlab("Relevant position") + 
      ylab("Effect size")+
      geom_hline(yintercept = 0, color = "black",linewidth=1) +
      theme_minimal()+
      theme(
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()
      )
  }
}

