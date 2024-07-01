import pandas as pd
import numpy as np
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import split_dimer



#-----------------------------------------------------
# helper dataset
#-----------------------------------------------------

# for match_cell_type_specificity()


df_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.hepG2.tab",sep="\t")
df_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.k562.tab",sep="\t")
df_dispersion_joint=pd.concat([df_dispersion_hepg2,df_dispersion_k562],axis=0).reset_index(drop=True)
df_dispersion_joint=df_dispersion_joint.drop_duplicates().reset_index(drop=True)

df_dispersion_hepg2["gini_rank"]=df_dispersion_hepg2["gini"].rank(ascending=False)
df_dispersion_k562["gini_rank"]=df_dispersion_k562["gini"].rank(ascending=False)

#-----------------------------------------------------
# helper functions 
#-----------------------------------------------------

color_mapping = {
    'sub': "#1f77b4", 
    'super': '#ff7f0e'
    }

def get_sub_super_tfs(suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_pca_coord_with_sub_super_{suffix}.csv")
    # select only TFs with very positive PC1
    # df=df[df["PC1"]>0.1].reset_index(drop=True)
    sub_tfs=df[df["tf_type"]=="sub"]["protein"].to_list()
    super_tfs=df[df["tf_type"]=="super"]["protein"].to_list()
    sub_tfs=split_dimer(sub_tfs)
    super_tfs=split_dimer(super_tfs)
    intersection=set(sub_tfs).intersection(set(super_tfs))
    sub_tfs=[tf for tf in sub_tfs if tf not in intersection]
    super_tfs=[tf for tf in super_tfs if tf not in intersection]
    return sub_tfs,super_tfs
        

def find_closest_numbers(array1, array2, num_closest=1):
    """
    For each item in array1, find the closest num_closest numbers in array2
    """
    array2_sorted = sorted(array2)
    all_closest_numbers = []
    for num1 in array1:
        closest_numbers = []
        for _ in range(num_closest):
            if not array2_sorted:
                break
            closest_num = min(array2_sorted, key=lambda x: abs(x - num1))
            closest_numbers.append(closest_num)
            array2_sorted.remove(closest_num)
        all_closest_numbers += closest_numbers
    return all_closest_numbers


def match_super_with_sub(sub_tfs,super_tfs,df_dispersion):
    sub_ranks=df_dispersion[df_dispersion["gene"].isin(sub_tfs)]["gini_rank"].to_list()
    super_ranks=df_dispersion[df_dispersion["gene"].isin(super_tfs)]["gini_rank"].to_list()
    matched_ranks=find_closest_numbers(sub_ranks,super_ranks,num_closest=1)
    matched_super_tfs=df_dispersion[df_dispersion["gini_rank"].isin(matched_ranks)]["gene"].to_list()
    return matched_super_tfs



def read_tfbs_maf(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_{file_suffix}.csv",header=None)
    df.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]
    df["dataset"]=file_suffix
    df["log10_AF"]=np.log10(df["AF"])
    return df



def read_maf_effect_size(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Effect_size_vs_maf/maf_with_effect_size_{file_suffix}.csv",header=None,index_col=0)
    df.reset_index(drop=True,inplace=True)
    df.columns=["Chromosome","Start","End","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
    return df



def join_tfbs_maf_effect_size(file_suffix):
    tfbs_maf=read_tfbs_maf(file_suffix)
    maf_effect_size=read_maf_effect_size(file_suffix)
    df=tfbs_maf.join(maf_effect_size.set_index(["Chromosome","Start","End","ID","REF","ALT","AF"]),on=["Chromosome","Start","End","ID","REF","ALT","AF"],how="inner")
    return df



#-----------------------------------------------------
# read in all SNPs, allele frequency, TFBS and effect size
#-----------------------------------------------------

# df_tfbs_maf_es is for effect size analysis
df_tfbs_maf_es=pd.concat([join_tfbs_maf_effect_size("enhancers_hepg2"),
                          join_tfbs_maf_effect_size("promoters_hepg2"),
                          join_tfbs_maf_effect_size("enhancers_k562"),
                          join_tfbs_maf_effect_size("promoters_k562")],
                          axis=0).reset_index(drop=True)
df_tfbs_maf_es=df_tfbs_maf_es[df_tfbs_maf_es["chip_evidence"]==True].reset_index(drop=True)
df_tfbs_maf_es["AF_bin"]=np.where(df_tfbs_maf_es["log10_AF"]<-4,"rare",np.where(df_tfbs_maf_es["log10_AF"]>-1,"common","low"))


df_tfbs_maf_es.drop(columns=["AF_bin"],inplace=True)

#--------------------------------------
# r variant enrichment per bin
#--------------------------------------

# TODO: summarize the enrichment sub in large effect size, in rare, low common variants

for cell_type in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562","hepg2","k562"]:
    if "hepg2" in cell_type:
        df_dispersion=df_dispersion_hepg2
        track_list = [0, 2, 4, 6]
    else:
        df_dispersion=df_dispersion_k562
        track_list = [1, 3, 5, 7]
    variant_type="rare"
    df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"].str.contains(cell_type)].reset_index(drop=True)
    df_subset = df_subset[df_subset["AF_bin"]==variant_type].reset_index(drop=True)
    sub_tfs,super_tfs=get_sub_super_tfs(cell_type)
    matched_super_tfs=match_super_with_sub(sub_tfs,super_tfs,df_dispersion) # no matching, no expected results.
    df_subset=df_subset[df_subset["motif"].isin(sub_tfs+matched_super_tfs)].reset_index(drop=True)
    df_subset["tf_type"]=np.where(df_subset["motif"].isin(sub_tfs),"sub","super")
    df_subset["tf_type"]=pd.Categorical(df_subset["tf_type"],categories=["sub","super"])
    for track_num in track_list:
        logger.info(f"Processing cell type {cell_type}, track {track_num}")
        df_subset_neg=df_subset[df_subset[f"track_{track_num}"]<0].reset_index(drop=True)
        df_subset_neg[f"track_{track_num}"]=df_subset_neg[f"track_{track_num}"].abs()
        sns.kdeplot(data=df_subset_neg,x=f"track_{track_num}",hue="tf_type",common_norm=False,cumulative=True)
        plt.title(f"cumulative distribution of effect size for {variant_type} variant\n{cell_type} track {track_num}")
        plt.axhline(y=1, color="black",linestyle=':')
        # add two vertical lines to show the 95 percentile effect size for sub and super TFs
        plt.axvline(x=df_subset_neg[df_subset_neg["tf_type"]=="sub"][f"track_{track_num}"].quantile(0.999),color=color_mapping["sub"],linestyle='--')
        plt.axvline(x=df_subset_neg[df_subset_neg["tf_type"]=="super"][f"track_{track_num}"].quantile(0.999),color=color_mapping["super"],linestyle='--')
        plt.savefig(f"Plots_{variant_type}_distribution/cumulative_distribution_{variant_type}_{cell_type}_track_{track_num}.pdf")
        plt.close()
        












#-----------------------------------------------------
# Archive
#-----------------------------------------------------

# bins=np.linspace(0,df_subset_neg[f"track_{track_num}"].max(),num=20)
# df_subset_neg=bin_and_label(df_subset_neg,f"track_{track_num}",bins)
# df_subset_neg["common"]=np.where(df_subset_neg["AF_bin"]=="common",1,0)
# df_subset_neg=df_subset_neg.groupby(["tf_type","Bin"]).agg({"common":"sum"}).reset_index()
# df_pivot=df_subset_neg.pivot_table(index="Bin",columns="tf_type",values="common",aggfunc="mean").reset_index()
# print(cell_type,track_num)
# print(df_pivot)







