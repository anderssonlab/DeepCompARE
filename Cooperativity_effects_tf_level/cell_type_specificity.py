import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from scipy.stats import pearsonr, spearmanr


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")

#-------------------------------
# calculate cell type specificity
#-------------------------------


def read_joint_dispersion():
    df_tf_dispersion_hepg2=pd.read_csv("TFs.dispersionEstimates.hepG2.tab",sep="\t")
    df_tf_dispersion_k562=pd.read_csv("TFs.dispersionEstimates.k562.tab",sep="\t")
    df_tf_dispersion=pd.concat([df_tf_dispersion_hepg2,df_tf_dispersion_k562],axis=0).reset_index(drop=True)
    df_tf_dispersion=df_tf_dispersion.drop_duplicates().reset_index(drop=True)
    df_tf_dispersion["gini_rank"]=df_tf_dispersion["gini"].rank(ascending=False)
    df_tf_dispersion["adjusted_dispersion_rank"]=df_tf_dispersion["adjusted_dispersion"].rank(ascending=False)
    df_tf_dispersion["mean_rank"]=df_tf_dispersion["mean"].rank(ascending=False)
    return df_tf_dispersion



tfs_redundant = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_redundant_lenient.txt",header=None)[0].to_list()
tfs_codependent = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_codependent_lenient.txt",header=None)[0].to_list()
df_tf_dispersion=read_joint_dispersion()
df_tf_dispersion["tf_type"]=np.where(df_tf_dispersion["gene"].isin(tfs_redundant),"redundant",np.where(df_tf_dispersion["gene"].isin(tfs_codependent),"codependent","other"))
df_tf_dispersion=df_tf_dispersion[df_tf_dispersion["tf_type"]!="other"].reset_index(drop=True)

tfs_codependent_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"]=="codependent"]["gini_rank"]
tfs_redundant_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"]=="redundant"]["gini_rank"]
stat, p_value = stats.ranksums(tfs_redundant_ranks,tfs_codependent_ranks)
print(f"p-value: {p_value}")
print(tfs_redundant_ranks.median()>tfs_codependent_ranks.median())

tfs_codependent_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"]=="codependent"]["mean_rank"]
tfs_redundant_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"]=="redundant"]["mean_rank"]
stat, p_value = stats.ranksums(tfs_redundant_ranks,tfs_codependent_ranks)
print(tfs_redundant_ranks.median()>tfs_codependent_ranks.median())

sns.kdeplot(data=df_tf_dispersion,x="gini",hue="tf_type",common_norm=False)
plt.title("Gini coefficient distribution")
plt.savefig("Plots/gini_distribution_lenient.pdf") 
plt.close()




#  calculate correlation between cooperativity ratio and gini

df_dispersion=read_joint_dispersion()
df_tf_cooperativity_ratio=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_cooperativity_ratio_lenient.csv")
df_dispersion=df_dispersion.merge(df_tf_cooperativity_ratio,left_on="gene",right_on="protein2",how="inner")
df_dispersion["tf_type"]=np.where(df_dispersion["protein2"].isin(tfs_redundant),"redundant",np.where(df_dispersion["protein2"].isin(tfs_codependent),"codependent","other"))

pcc,p_pear=pearsonr(df_dispersion["cooperativity_ratio"],df_dispersion["gini"])
scc,p_spear=spearmanr(df_dispersion["cooperativity_ratio"],df_dispersion["gini"])

# sns scatter plot, x=cooperativity_ratio, y=gini,shape by tf_type, add a linear fit
sns.regplot(data=df_dispersion,x="cooperativity_ratio",y="gini",scatter=False,color="black")
sns.scatterplot(data=df_dispersion,x="cooperativity_ratio",y="gini",hue="tf_type")
plt.ylim(0,1.2)
# annotate with spearman and pearson correlation, retain 2 decimal points
plt.text(0.5, 1.1, f"Pearson corr: {pcc:.2f} (p: {p_pear:.2f})\nSpearman corr: {scc:.2f} (p: {p_spear:.2f})")
plt.title("Cooperativity ratio vs Gini coefficient")
plt.savefig("Plots/cooperativity_ratio_vs_gini_lenient.pdf")
plt.close()

