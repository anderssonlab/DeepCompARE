import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from loguru import logger
from scipy import stats

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
sub_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/sub_tfs.txt",header=None)[0].to_list()
super_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/super_tfs.txt",header=None)[0].to_list()
df_tf_individual_effect=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect.csv")
df_tf_individual_effect["tf_type"]=np.where(df_tf_individual_effect["protein"].isin(sub_tfs),"sub",np.where(df_tf_individual_effect["protein"].isin(super_tfs),"super","other"))



df_tf_individual_effect[df_tf_individual_effect["tf_type"]=="sub"]
df_tf_individual_effect[df_tf_individual_effect["tf_type"]=="super"]