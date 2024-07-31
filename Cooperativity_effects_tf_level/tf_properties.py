import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from scipy.stats import pearsonr, spearmanr


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")

#-------------------------------
# calculate cell type specificity
#-------------------------------


def read_joint_dispersion():
    df_tf_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/TFs.dispersionEstimates.hepG2.tab",sep="\t")
    df_tf_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/TFs.dispersionEstimates.k562.tab",sep="\t")
    df_tf_dispersion=pd.concat([df_tf_dispersion_hepg2,df_tf_dispersion_k562],axis=0).reset_index(drop=True)
    df_tf_dispersion=df_tf_dispersion.drop_duplicates().reset_index(drop=True)
    df_tf_dispersion["gini_rank"]=df_tf_dispersion["gini"].rank(ascending=False)
    df_tf_dispersion["adjusted_dispersion_rank"]=df_tf_dispersion["adjusted_dispersion"].rank(ascending=False)
    df_tf_dispersion["mean_rank"]=df_tf_dispersion["mean"].rank(ascending=False)
    return df_tf_dispersion



tfs_redundant = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tfs_redundant.txt",header=None)[0].to_list()
tfs_codependent = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tfs_codependent.txt",header=None)[0].to_list()
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
plt.savefig("Plots/gini_distribution.pdf") 
plt.close()




#  calculate correlation between cooperativity ratio and gini

df_dispersion=read_joint_dispersion()
df_tf_cooperativity_ratio=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_cooperativity_ratio.csv")
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
plt.savefig("Plots/cooperativity_ratio_vs_gini.pdf")
plt.close()


#-------------------------------
# plot individual TF effect v.s. cr
#-------------------------------


df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_ism_dstat_by_track_chip_true.csv")
df_coop=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_cooperativity_ratio.csv")
# merge
df=df.merge(df_coop,left_on="protein",right_on="protein2",how="inner")
# sns scatter plot
#for track in ["cage","sure"]:
for track in ["cage","dhs","starr","sure"]:
    sns.regplot(data=df,x=f"{track}",y="cooperativity_ratio")
    plt.title(f"{track} vs cooperativity ratio")
    # annotate with correlation
    pcc,p_pear=pearsonr(df[track],df["cooperativity_ratio"])
    scc,p_spear=spearmanr(df[track],df["cooperativity_ratio"])
    plt.text(0.3, 0.9, f"Pearson R: {pcc:.2f} (p: {p_pear:.2f})\nSpearman R: {scc:.2f} (p: {p_spear:.2f})")
    # annotate with protein name
    for i in range(len(df)):
        plt.text(df[track][i],df["cooperativity_ratio"][i],df["protein"][i],fontsize=6)
    plt.savefig(f"Plots/cr_vs_{track}_effect_dstat.pdf")
    plt.close()


#-------------------------------
# Color subs and tfs_codependents in PCA
#-------------------------------

# TODO: refactor code to use /isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_pca_coord_with_sub_tfs_codependent_hepg2.csv
def plot_pca_color_by_sub_tfs_codependent(df,title,outname):
    features=["cage","dhs","starr","sure"]
    features_prefixed=[f"dstat_ism_{x}" for x in features]
    x = df.loc[:,features_prefixed].values
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    pca_components = pca.components_
    pc1_label = " ".join(f"{coef:.2f}{feat}" for coef, feat in zip(pca_components[0], features))
    pc2_label = " ".join(f"{coef:.2f}{feat}" for coef, feat in zip(pca_components[1], features))
    explained_variance = pca.explained_variance_ratio_
    pc1_explained = f"{explained_variance[0]:.2%}"
    pc2_explained = f"{explained_variance[1]:.2%}"
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    principalDf["tf_type"]=df["tf_type"]
    principalDf["alpha_value"]=np.where(principalDf["tf_type"].isin(["sub","tfs_codependent"]),1,0.1)
    plt.figure(figsize=(6, 6))
    sns.scatterplot(x='PC1', y='PC2', data=principalDf, hue='tf_type', alpha=principalDf["alpha_value"], s=10)
    plt.xlabel(f"PC1 = {pc1_label} (EV: {pc1_explained})")
    plt.ylabel(f"PC2 = {pc2_label} (EV: {pc2_explained})")
    plt.title(title)
    plt.savefig(outname)
    plt.close()
    
    
def plot_pca_color_by_sub_tfs_codependent_all(cell_type):
    df=read_tf_individual_effect()
    df=df[df["dataset"]==cell_type].reset_index(drop=True)
    tfs_redundant,tfs_codependent_tfs=get_sub_tfs_codependent_annotation(cell_type)
    df["tf_type"]=np.where(df["protein"].isin(tfs_redundant),"sub",np.where(df["protein"].isin(tfs_codependent_tfs),"tfs_codependent","other"))
    plot_pca_color_by_sub_tfs_codependent(df,cell_type,f"pca_{cell_type}_color_sub_tfs_codependent.pdf")
    
    
plot_pca_color_by_sub_tfs_codependent_all("hepg2")
plot_pca_color_by_sub_tfs_codependent_all("k562")


df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_pca_coord_k562_chip_true.csv")
# add cooperativity ratio
df_cooperativity_ratio=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_cooperativity_ratio.csv")
# merge
df=df.merge(df_cooperativity_ratio,left_on="protein",right_on="protein2",how="inner")
# sns scatter plot, x=PC1, y=PC2, hue=cooperativity_ratio
sns.scatterplot(data=df,x="PC1",y="PC2",hue="cooperativity_ratio")
plt.title("PC1 vs PC2, color by cooperativity ratio")
# annotate with protein name
for i in range(len(df)):
    plt.text(df["PC1"][i],df["PC2"][i],df["protein"][i],fontsize=4)

plt.savefig("Plots/pca_k562_chip_true_color_cooperativity_ratio.pdf")
plt.close()