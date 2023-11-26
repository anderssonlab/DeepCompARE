

df_plot <- data.frame(MPRA_effect_size=df_truth$Value,
                      Delta_z=df_pred$z_cage_hepg2_alt-df_pred$z_cage_hepg2_ref)

df_plot <- data.frame(MPRA_effect_size=df_truth$Value,
                      Delta_z=log2(exp(df_pred$log1p_signal_cage_hepg2_alt-1)/exp(df_pred$log1p_signal_cage_hepg2_ref-1)))
cor(df_plot$MPRA_effect_size,df_plot$Delta_z)
ggplot(df_plot,aes(x=Delta_z,y=MPRA_effect_size))+
  geom_point()+
  xlab("Delta_z_by_classification")+
  geom_smooth(method="lm", se=FALSE, col="blue")+
  annotate(geom="text", x=-0.5, y=1, label="PCC=0.755", size=4.5, col="black") +
  theme_bw()


df_truth <- read.csv("Validation_SNP_MPRA/MPRA_Promoter_HepG2_F9.csv")
df_pred <- read.csv("Validation_SNP_MPRA/AstigCRConv5D_Dataset_final_CR_Trainer/log1p_signal_Promoter_HepG2_F9.csv")
df_plot <- data.frame(MPRA_effect_size=df_truth$Value,
                      Delta_z=df_pred$log1p_signal_cage_hepg2_alt-df_pred$log1p_signal_cage_hepg2_ref)

df_plot <- data.frame(MPRA_effect_size=df_truth$Value,
                      Delta_z=log2(exp(df_pred$log1p_signal_cage_hepg2_alt-1)/exp(df_pred$log1p_signal_cage_hepg2_ref-1)))
cor(df_plot$MPRA_effect_size,df_plot$Delta_z)
ggplot(df_plot,aes(x=Delta_z,y=MPRA_effect_size))+
  geom_point()+
  xlab("Delta_signal_by_regression")+
  geom_smooth(method="lm", se=FALSE, col="blue")+
  annotate(geom="text", x=0, y=1, label="PCC=0.414", size=4.5, col="black") +
  theme_bw()
