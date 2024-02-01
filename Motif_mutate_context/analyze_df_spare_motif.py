import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

import sys
sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from plotting import scatter_plot_with_annotation

df=pd.read_csv("df_spare_motif_cage_k562.csv")
df.columns
df.shape

#-------------------------
# Step 1: clean up weired values
#-------------------------
# # 1. negative feat_imp_orig
# df_neg_feat_imp=df[df["feat_imp_orig"]<0] # 2276
# sns.scatterplot(data=df_neg_feat_imp,x="feat_imp_orig",y="feat_imp_mut")
# plt.title("Negative feat_imp_orig")
# plt.savefig("Negative_feat_imp_orig.png")
# plt.close()

# # 2. negative feat_imp_mut
# df_neg_feat_imp=df[df["feat_imp_mut"]<0] # 46, seems like noise, not pattern
# sns.scatterplot(data=df_neg_feat_imp,x="feat_imp_orig",y="feat_imp_mut")
# plt.title("Negative feat_imp_mut")
# plt.savefig("Negative_feat_imp_mut.png")
# plt.close()

df=df[(df["feat_imp_orig"]>=0) & (df["feat_imp_mut"]>=0)] #3448081


#-------------------------
# Step 2: add info
#-------------------------
# add TF family information and calculate percentage remain
df_tf_family=pd.read_csv("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2024_CORE_extracted_tfs.csv",
                         usecols=["ID","tf_family"])
df_tf_family["ID"]=df_tf_family["ID"].str.upper()


df["protein"].isin(df_tf_family["ID"]).sum() # 3448081
df.shape # 3448081
df=df.merge(df_tf_family,left_on="protein",right_on="ID")
df=df.drop(columns=["ID"])

df["percentage_remain"]=df["feat_imp_mut"]/df["feat_imp_orig"]



#--------------------------------
# Step 3: protein level analysis
#--------------------------------

# group by TF, for each TF, compute the median percentage_remain, and mean feat_imp_orig  and mean feat_imp_mut
df_protein_percentage_remain=df.groupby("protein").agg({"percentage_remain":"median","feat_imp_orig":"mean","feat_imp_mut":"mean"}).reset_index()

# Spearman correlation with ChIP-Seq KS test
df_protein_percentage_remain.sort_values(by="protein",ascending=False,inplace=True)
df_ks=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/ks_stat_featimp_and_motif_score.csv")
df_ks.sort_values(by="protein",ascending=False,inplace=True)
assert np.all(df_ks.protein==df_protein_percentage_remain.protein)
df_protein_percentage_remain["rank"]=df_protein_percentage_remain["percentage_remain"].rank(method="dense",ascending=False)
df_ks["rank"]=df_ks["feat_imp_d_stat"].rank(method="dense",ascending=False)
corr, p_value = spearmanr(df_protein_percentage_remain["rank"], df_ks["rank"]) # -0.24, 0.0004



scatter_plot_with_annotation(df_protein_percentage_remain, 
                             "feat_imp_mut", 
                             "percentage_remain",
                             "protein", 
                             "Plots_spare_motif/protein_level_featimp_mut_vs_percentage_remain.png",
                             xlab="mean_feat_imp_mut", ylab="median_percentage_remain")

scatter_plot_with_annotation(df_protein_percentage_remain,
                             "feat_imp_orig", 
                             "percentage_remain",
                             "protein", 
                             "Plots_spare_motif/protein_level_featimp_orig_vs_percentage_remain.png",
                             xlab="mean_feat_imp_orig", ylab="median_percentage_remain")


scatter_plot_with_annotation(df_protein_percentage_remain,
                             "feat_imp_orig", 
                             "feat_imp_mut",
                             "protein", 
                             "Plots_spare_motif/protein_level_featimp_orig_vs_featimp_mut.png",
                             xlab="mean_feat_imp_orig", ylab="mean_feat_imp_mut")


#--------------------------------
# Step 4: TF family level analysis
#--------------------------------
df_tf_family_percentage_remain=df.groupby("tf_family").agg({"percentage_remain":"median","feat_imp_orig":"mean","feat_imp_mut":"mean"}).reset_index()

scatter_plot_with_annotation(df_tf_family_percentage_remain,
                             "feat_imp_orig", 
                             "percentage_remain",
                             "tf_family", 
                             "Plots_spare_motif/tf_family_level_featimp_orig_vs_percentage_remain.png",
                             xlab="mean_feat_imp_orig", ylab="median_percentage_remain")
scatter_plot_with_annotation(df_tf_family_percentage_remain, 
                             "feat_imp_mut", 
                             "percentage_remain",
                             "tf_family", 
                             "Plots_spare_motif/tf_family_level_featimp_mut_vs_percentage_remain.png",
                             xlab="mean_feat_imp_mut", ylab="median_percentage_remain")
scatter_plot_with_annotation(df_tf_family_percentage_remain,
                             "feat_imp_orig", 
                             "feat_imp_mut",
                             "tf_family", 
                             "Plots_spare_motif/tf_family_level_featimp_orig_vs_featimp_mut.png",
                             xlab="mean_feat_imp_orig", ylab="mean_feat_imp_mut")



#--------------------------------
# Step 5: ETS family analysis 
#--------------------------------
# select tf_family containing "Ets-related"
df_ets=df[df["tf_family"].str.contains("Ets-related").fillna(False)]
df_ets_percentage_remain=df_ets.groupby("protein").agg({"percentage_remain":"median","feat_imp_orig":"mean","feat_imp_mut":"mean"}).reset_index()
scatter_plot_with_annotation(df_ets_percentage_remain,
                                "feat_imp_orig", 
                                "percentage_remain",
                                "protein", 
                                "Plots_spare_motif/ETS_family_level_featimp_orig_vs_percentage_remain.png",
                                xlab="mean_feat_imp_orig", ylab="median_percentage_remain")

scatter_plot_with_annotation(df_ets_percentage_remain,
                                "feat_imp_mut", 
                                "percentage_remain",
                                "protein", 
                                "Plots_spare_motif/ETS_family_level_featimp_mut_vs_percentage_remain.png",
                                xlab="mean_feat_imp_mut", ylab="median_percentage_remain")

scatter_plot_with_annotation(df_ets_percentage_remain,
                                "feat_imp_orig", 
                                "feat_imp_mut",
                                "protein", 
                                "Plots_spare_motif/ETS_family_level_featimp_orig_vs_featimp_mut.png",
                                xlab="mean_feat_imp_orig", ylab="mean_feat_imp_mut")