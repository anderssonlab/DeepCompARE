import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

import sys
sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from plotting import scatter_plot_with_annotation

df_promoters_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_promoters_k562.csv")
df_promoters_k562["file_name"]="promoters_k562"

df_promoters_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_promoters_hepg2.csv")
df_promoters_hepg2["file_name"]="promoters_hepg2"

df_enhancers_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_enhancers_k562.csv")
df_enhancers_k562["file_name"]="enhancers_k562"

df=pd.concat([df_promoters_k562,df_promoters_hepg2,df_enhancers_k562],axis=0)


#-------------------------
# Step 1: clean up weired values
#-------------------------
# 1. negative feat_imp_orig
df_neg_feat_imp=df[df["feat_imp_orig"]<0] # 770
# count occurances of each protein
df_neg_feat_imp["protein"].value_counts() # No.1 ZNF384 378
df["protein"].value_counts() # ZNF384 13993

# 2. negative feat_imp_mut
df_neg_feat_imp_remove_context=df[df["feat_imp_remove_context"]<0] # 46 
df_neg_feat_imp_remove_context["protein"].value_counts() # No.1 ZNF384 23



df=df[(df["feat_imp_orig"]>=0) & (df["feat_imp_remove_context"]>=0)] 
df["log_ratio"]=np.log2(df["feat_imp_remove_context"]/df["feat_imp_orig"])


#-------------------------
# Step 2: add tf families
#-------------------------
# add TF family information and calculate percentage remain
# df_tf_family=pd.read_csv("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2024_CORE_extracted_tfs.csv",
#                          usecols=["ID","tf_family"])
# df_tf_family["ID"]=df_tf_family["ID"].str.upper()


# df["protein"].isin(df_tf_family["ID"]).sum() # 3448081
# df.shape # 3448081
# df=df.merge(df_tf_family,left_on="protein",right_on="ID")
# df=df.drop(columns=["ID"])



#--------------------------------
# Step 3: protein level analysis
#--------------------------------

# group by TF, for each TF, compute the median percentage_remain, and mean feat_imp_orig  and mean feat_imp_mut
df_log_ratio=df.groupby("protein").agg({"log_ratio":"mean","feat_imp_orig":"mean","feat_imp_remove_context":"mean"}).reset_index()



df_ks_promoters_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Motif1_highlight_chip/summary_and_ks_test_promoters_k562.csv")
df_ks_promoters_k562["file_name"]="promoters_k562"

df_ks_promoters_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Motif1_highlight_chip/summary_and_ks_test_promoters_hepg2.csv")
df_ks_promoters_hepg2["file_name"]="promoters_hepg2"

df_ks_enhancers_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Motif1_highlight_chip/summary_and_ks_test_enhancers_k562.csv")
df_ks_enhancers_k562["file_name"]="enhancers_k562"

df_ks=pd.concat([df_ks_promoters_k562,df_ks_promoters_hepg2,df_ks_enhancers_k562],axis=0)
df_ks.reset_index(drop=True,inplace=True)


# join df_log_ratio and df_ks by protein




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


