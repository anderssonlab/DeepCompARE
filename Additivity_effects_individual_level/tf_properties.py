import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from loguru import logger
from scipy import stats
from sklearn.decomposition import PCA
from adjustText import adjust_text
#-------------------------------
# calculate cell type specificity
#-------------------------------

def sub_sup_gini(rank="gini_rank"):
    sub_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/sub_tfs.txt",header=None)[0].to_list()
    print(sub_tfs)
    super_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/super_tfs.txt",header=None)[0].to_list()
    print(super_tfs)
    df_tf_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.hepG2.tab",sep="\t")
    df_tf_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.k562.tab",sep="\t")
    df_tf_dispersion=pd.concat([df_tf_dispersion_hepg2,df_tf_dispersion_k562],axis=0).reset_index(drop=True)
    df_tf_dispersion=df_tf_dispersion.drop_duplicates().reset_index(drop=True)
    # get ranks
    df_tf_dispersion["gini_rank"]=df_tf_dispersion["gini"].rank(ascending=False)
    df_tf_dispersion["adjusted_dispersion_rank"]=df_tf_dispersion["adjusted_dispersion"].rank(ascending=False)
    # add tf type: sub or super
    df_tf_dispersion["tf_type"]=np.where(df_tf_dispersion["gene"].isin(sub_tfs),"sub",np.where(df_tf_dispersion["gene"].isin(super_tfs),"super","other"))
    # cell type specificity using wilcoxon rank sum test
    # for subs
    sub_tfs_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"]=="sub"][rank]
    other_tfs_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"].isin(["super","other"])][rank]
    stat, p_value = stats.ranksums(sub_tfs_ranks, other_tfs_ranks) # 0.03
    print("sub",p_value)
    print(sub_tfs_ranks.median()>other_tfs_ranks.median()) # True
    # for super
    super_tfs_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"]=="super"][rank]
    other_tfs_ranks = df_tf_dispersion[df_tf_dispersion["tf_type"].isin(["sub","other"])][rank]
    stat, p_value = stats.ranksums(super_tfs_ranks, other_tfs_ranks)  # 0.51
    print("super",p_value)
    print(super_tfs_ranks.median()<other_tfs_ranks.median()) # True
    stat, p_value = stats.ranksums(sub_tfs_ranks, super_tfs_ranks)  # 0.03
    print("sub vs super",p_value)
    # plot the distribution of gini
    sns.kdeplot(data=df_tf_dispersion,x="gini",hue="tf_type",common_norm=False)
    plt.savefig("Plots/gini_tf_type_distribution.pdf") 
    plt.close()
    
    
    
sub_sup_gini()

#-------------------------------
# plot individual TF effect
#-------------------------------
def read_tf_individual_effect():
    df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_cell_type.csv")
    columns=[col for col in df.columns if "dstat_ism" in col]
    df=df[columns+["protein","dataset"]]
    return df

def get_sub_super_annotation(cell_type):
    df_tf_property=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_additivity_property.csv",index_col=0)
    datasets=[f"enhancers_{cell_type}", f"promoters_{cell_type}"]
    df_tf_property=df_tf_property[datasets]
    df_tf_property["sub_counts"]=df_tf_property.apply(lambda x: x.value_counts().get("sub",0),axis=1)
    df_tf_property["super_counts"]=df_tf_property.apply(lambda x: x.value_counts().get("super",0),axis=1)
    sub_tfs =  df_tf_property[(df_tf_property["sub_counts"]>0) & (df_tf_property["super_counts"]==0)].index.to_list()
    super_tfs = df_tf_property[(df_tf_property["super_counts"]>0) & (df_tf_property["sub_counts"]==0)].index.to_list()
    return sub_tfs,super_tfs

def plot_sub_super_tf_effect_by_track(cell_type):
    sub_tfs,super_tfs=get_sub_super_annotation(cell_type)
    # read individual effect
    df=read_tf_individual_effect()
    df=df[df["dataset"]==cell_type].reset_index(drop=True)
    df["tf_type"]=np.where(df["protein"].isin(sub_tfs),"sub",np.where(df["protein"].isin(super_tfs),"super","other"))
    # select only rows with sub or super 
    df=df[df["tf_type"].isin(["sub","super"])].reset_index(drop=True)
    # pivot df to long format
    columns=[col for col in df.columns if "dstat_ism" in col]
    df=pd.melt(df,id_vars=["protein","tf_type"],value_vars=columns,var_name="track",value_name="dstat_ism")
    # remove "dstat_ism_" from track
    df["track"]=df["track"].str.replace("dstat_ism_","")
    sns.violinplot(data=df,x="track",y="dstat_ism",hue="tf_type",split=True,inner="quartile")
    plt.title(f"{cell_type}")
    plt.savefig(f"sub_super_tf_effect_{cell_type}.pdf")
    plt.close()


plot_sub_super_tf_effect_by_track("hepg2")
plot_sub_super_tf_effect_by_track("k562")



#-------------------------------
# Color subs and supers in PCA
#-------------------------------

# TODO: refactor code to use /isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_pca_coord_with_sub_super_hepg2.csv
def plot_pca_color_by_sub_super(df,title,outname):
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
    principalDf["alpha_value"]=np.where(principalDf["tf_type"].isin(["sub","super"]),1,0.1)
    plt.figure(figsize=(6, 6))
    sns.scatterplot(x='PC1', y='PC2', data=principalDf, hue='tf_type', alpha=principalDf["alpha_value"], s=10)
    plt.xlabel(f"PC1 = {pc1_label} (EV: {pc1_explained})")
    plt.ylabel(f"PC2 = {pc2_label} (EV: {pc2_explained})")
    plt.title(title)
    plt.savefig(outname)
    plt.close()
    
    
def plot_pca_color_by_sub_super_all(cell_type):
    df=read_tf_individual_effect()
    df=df[df["dataset"]==cell_type].reset_index(drop=True)
    sub_tfs,super_tfs=get_sub_super_annotation(cell_type)
    df["tf_type"]=np.where(df["protein"].isin(sub_tfs),"sub",np.where(df["protein"].isin(super_tfs),"super","other"))
    plot_pca_color_by_sub_super(df,cell_type,f"pca_{cell_type}_color_sub_super.pdf")
    
    
plot_pca_color_by_sub_super_all("hepg2")
plot_pca_color_by_sub_super_all("k562")


