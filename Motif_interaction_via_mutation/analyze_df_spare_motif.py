import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from scipy.stats import spearmanr

df=pd.read_csv("df_spare_motif_cage_k562.csv")
df.columns
df.shape

# clean up weired values
# 1. negative feat_imp_orig
df_neg_feat_imp=df[df["feat_imp_orig"]<0] # 2276
sns.scatterplot(data=df_neg_feat_imp,x="feat_imp_orig",y="feat_imp_mut")
plt.title("Negative feat_imp_orig")
plt.savefig("Negative_feat_imp_orig.png")
plt.close()
# 2. negative feat_imp_mut
df_neg_feat_imp=df[df["feat_imp_mut"]<0] # 46, seems like noise, not patter
sns.scatterplot(data=df_neg_feat_imp,x="feat_imp_orig",y="feat_imp_mut")
plt.title("Negative feat_imp_mut")
plt.savefig("Negative_feat_imp_mut.png")
plt.close()


df=df[(df["feat_imp_orig"]>=0) & (df["feat_imp_mut"]>=0)] #3448081
df["percentage_remain"]=df["feat_imp_mut"]/df["feat_imp_orig"]
# group by TF, for each TF, compute the median percentage_remain, and mean feat_imp_orig  and mean feat_imp_mut
df_median_percentage_remain=df.groupby("protein").agg({"percentage_remain":"median","feat_imp_orig":"mean","feat_imp_mut":"mean"}).reset_index()
df_median_percentage_remain.sort_values(by="protein",ascending=False,inplace=True)


df_ks=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/ks_stat_featimp_and_motif_score.csv")
df_ks.sort_values(by="protein",ascending=False,inplace=True)
assert np.all(df_ks.protein==df_median_percentage_remain.protein)
df_median_percentage_remain["rank"]=df_median_percentage_remain["percentage_remain"].rank(method="dense",ascending=False)
df_ks["rank"]=df_ks["feat_imp_d_stat"].rank(method="dense",ascending=False)
corr, p_value = spearmanr(df_median_percentage_remain["rank"], df_ks["rank"]) # -0.24, 0.0004




df_median_percentage_remain.sort_values(by="percentage_remain",ascending=True,inplace=True)
df_ks.sort_values(by="feat_imp_d_stat",ascending=False,inplace=True)
df_ks.iloc[0:10,[0,1,4,5,6]]

texts = []
plt.figure(figsize=(10, 8))  
scatter_plot = sns.scatterplot(data=df_median_percentage_remain, x="feat_imp_orig", y="percentage_remain")
for line in range(0, df_median_percentage_remain.shape[0]):
    text=scatter_plot.text(df_median_percentage_remain.feat_imp_orig[line], df_median_percentage_remain.percentage_remain[line], 
                            df_median_percentage_remain.protein[line], horizontalalignment='left', 
                            size='small', color='black')
    texts.append(text)
adjust_text(texts)
plt.xlabel("mean_feat_imp_orig")
plt.ylabel("median_percentage_remain")
plt.savefig("feat_imp_orig_vs_percentage_remain.png")
plt.close()



texts = []
plt.figure(figsize=(10, 8))  
scatter_plot = sns.scatterplot(data=df_median_percentage_remain, x="feat_imp_mut", y="percentage_remain")
for line in range(0, df_median_percentage_remain.shape[0]):
    text=scatter_plot.text(df_median_percentage_remain.feat_imp_mut[line], df_median_percentage_remain.percentage_remain[line], 
                            df_median_percentage_remain.protein[line], horizontalalignment='left', 
                            size='small', color='black')
    texts.append(text)
adjust_text(texts)
plt.xlabel("mean_feat_imp_mut")
plt.ylabel("median_percentage_remain")
plt.savefig("feat_imp_mut_vs_percentage_remain.png")
plt.close()