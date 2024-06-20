import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from loguru import logger
from scipy import stats
import random


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_and_label, calc_or, plot_or



#-----------------------------------------------------
# helper dataset
#-----------------------------------------------------

# for match_cell_type_specificity()


df_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.hepG2.tab",sep="\t")
df_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.k562.tab",sep="\t")
df_dispersion_joint=pd.concat([df_dispersion_hepg2,df_dispersion_k562],axis=0).reset_index(drop=True)
df_dispersion_joint=df_dispersion_joint.drop_duplicates().reset_index(drop=True)




#-----------------------------------------------------
# helper functions 
#-----------------------------------------------------

color_mapping = {
    'sub': "#1f77b4", 
    'super': '#ff7f0e'
    }

def get_contextual_tf_type(df,dataset):
    df_tf_property=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_additivity_property.csv")
    sub_tfs = df_tf_property[df_tf_property[dataset]=="sub"]["protein"].to_list()
    super_tfs = df_tf_property[df_tf_property[dataset]=="super"]["protein"].to_list()
    tf_type= np.where(df["motif"].isin(sub_tfs),"sub",np.where(df["motif"].isin(super_tfs),"super","other"))
    # tf_type=pd.Categorical(tf_type,categories=["sub","super","other"])
    return tf_type


def generate_random_shift(n):
    return [random.choice([1, -1]) for _ in range(n)]


def find_closest_numbers(array1, array2, num_closest):
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


def match_cell_type_specificity(df,df_dispersion,tf_type,num_closest=1):
    df_dispersion["gini_rank"]=df_dispersion["gini"].rank(ascending=False)
    target_tfs = df[df["TF_type"] == tf_type]["motif"].unique().tolist()
    other_tfs = df[df["TF_type"] == "other"]["motif"].unique().tolist()
    if len(target_tfs)>=len(other_tfs):
        return df[df["motif"].isin(target_tfs+other_tfs)].reset_index(drop=True).copy()
    target_ranks=df_dispersion[df_dispersion["gene"].isin(target_tfs)]["gini_rank"].to_list()
    other_ranks=df_dispersion[df_dispersion["gene"].isin(other_tfs)]["gini_rank"].to_list()
    matched_ranks=find_closest_numbers(target_ranks,other_ranks,num_closest)
    matched_tfs=df_dispersion[df_dispersion["gini_rank"].isin(matched_ranks)]["gene"].to_list()
    df_subset=df[df["motif"].isin(target_tfs+matched_tfs)].reset_index(drop=True).copy()
    return df_subset


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
df_tfbs_maf_es["AF_bin"]=np.where(df_tfbs_maf_es["log10_AF"]<-3,"rare",np.where(df_tfbs_maf_es["log10_AF"]<-2,"low","common"))
df_tfbs_maf_es=df_tfbs_maf_es[df_tfbs_maf_es["AF_bin"]!="low"].reset_index(drop=True)


sub_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/sub_tfs.txt",header=None)[0].to_list()
super_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/super_tfs.txt",header=None)[0].to_list()


# df_tfbs_maf_es.drop(columns=["TF_type"],inplace=True)
#--------------------------------------------------------------------------------------------
# 1. Allele frequency analysis, sub v.s. rest
#--------------------------------------------------------------------------------------------

df_tfbs_maf_es["TF_type"]=np.where(df_tfbs_maf_es["motif"].isin(sub_tfs),"sub","other")

# plot distribution of log10_AF, hue by TF_type
sns.kdeplot(data=df_tfbs_maf_es,x="log10_AF",hue="TF_type",common_norm=False)
plt.title("Allele frequency distribution of SNVs overlapping REs")
plt.savefig("Plots/af_distribution.pdf")
plt.close()

sns.histplot(data=df_tfbs_maf_es,x="log10_AF",hue="TF_type",bins=100)
plt.title(f"Allele frequency of SNVs overlapping all REs")
plt.yscale("log")
plt.savefig(f"Plots/af_distribution_grouped_all.pdf")
plt.close()

sns.histplot(data=df_tfbs_maf_es,x="log10_AF",hue="TF_type",bins=100,multiple="fill")
plt.title(f"Allele frequency of SNVs overlapping all REs")
plt.savefig(f"Plots/af_distribution_grouped_fill_all.pdf")
plt.close()


df_binned = df_tfbs_maf_es.loc[:, ["AF_bin", "TF_type"]].copy().groupby(["AF_bin", "TF_type"]).size().unstack()
df_plot= calc_or(df_binned, "AF_bin", "TF_type")
plot_or(df_plot,"AF_bin", "odds_ratio", "TF_type",
        "Enrichment analysis, sub v.s. rest",{'sub': "#1f77b4"},
        out_name="Plots/af_enrichment_not_matched.pdf",
        ax=None) # Comparing to rest, subs are depleted in common

# match by cell type specificity
df_matched=match_cell_type_specificity(df_tfbs_maf_es,df_dispersion_joint,"sub",num_closest=3)
df_matched_binned = df_matched.loc[:, ["AF_bin", "TF_type"]].copy().groupby(["AF_bin", "TF_type"]).size().unstack()
df_matched_plot= calc_or(df_matched_binned, "AF_bin", "TF_type")
plot_or(df_matched_plot,"AF_bin", "odds_ratio", "TF_type",
        "Enrichment analysis, sub v.s. rest",{'sub': "#1f77b4"},
        out_name="Plots/af_enrichment_matched.pdf",
        ax=None) # After cell type specificity matching, subs are less depleted in common

#--------------------------------------------------------------------------------------------
# 2. Stratify by AF,match cell type specificity plot enrichment of sub/super over effect size bins
#--------------------------------------------------------------------------------------------

def get_or_for_plot(df,x_colname,hue_colname="TF_type"):
    df = df.loc[:, [x_colname, hue_colname]].copy().groupby([x_colname, hue_colname]).size().unstack()
    df = df.fillna(0)
    df_or = calc_or(df, x_colname, hue_colname, out_group="other")
    return df_or


def plot_stratified_enrichment(df,df_dispersion,strata_colname,strata_prefix, x_colname, y_colname, sup_title, out_name, color_mapping, match):
    strata_bins = df[strata_colname].unique()
    fig, axs = plt.subplots(len(strata_bins), 1, figsize=(8, 8), sharex=True)
    for i, strata_bin in enumerate(strata_bins):
        # logger.info(f"Processing {strata_prefix}: {strata_bin}")
        df_strata_subset = df[df[strata_colname] == strata_bin].reset_index(drop=True)
        if df_strata_subset.shape[0] < 10:
            continue
        if match:
            df_strata_subset_super_matched=match_cell_type_specificity(df_strata_subset,df_dispersion,"super",num_closest=1)
            df_strata_subset_sub_matched=match_cell_type_specificity(df_strata_subset,df_dispersion,"sub",num_closest=3)
            df_plot_sub = get_or_for_plot(df_strata_subset_sub_matched, x_colname)
            df_plot_super = get_or_for_plot(df_strata_subset_super_matched, x_colname)
            df_plot=pd.concat([df_plot_sub,df_plot_super],axis=0).reset_index(drop=True)
        else:
            df_strata_subset = df_strata_subset[df_strata_subset["TF_type"] != "other"].reset_index(drop=True)
            df_plot=get_or_for_plot(df_strata_subset, x_colname)
        plot_or(df_plot, x_colname, y_colname, "TF_type",
                f"{strata_prefix}: {strata_bin}", color_mapping, ax=axs[i])
    axs[-1].set_xlabel(x_colname)
    fig.text(0.04, 0.5, y_colname, va='center', rotation='vertical')
    plt.suptitle(sup_title) 
    plt.tight_layout()
    plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.15)
    plt.savefig(out_name)
    plt.close()
    



for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
    df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
    df_subset["TF_type"]=get_contextual_tf_type(df_subset,dataset)
    if "hepg2" in dataset:
        track_list = [0, 2, 4, 6]
        df_dispersion=df_dispersion_hepg2
    else:
        track_list = [1, 3, 5, 7]
        df_dispersion=df_dispersion_k562
    for track_num in track_list:
        df_subset = bin_and_label(df_subset, f"track_{track_num}", [-np.inf, -0.1, 0.1, np.inf],"es_bin")
        for match in [False,True]:
            logger.info(f"Processing dataset {dataset}, track {track_num}, match {match}")
        # stratify by AF_bin, plot es_enrichment
            plot_stratified_enrichment(df_subset,df_dispersion,
                                    "AF_bin","Variant type","es_bin", "odds_ratio", 
                                    f"{dataset}, track {track_num}", 
                                    f"Plots_enrichment/match_{match}_es_enrichment_{dataset}_track{track_num}.pdf",
                                    color_mapping,
                                    match=match)
            # stratify by es_bin, plot af_enrichment
            plot_stratified_enrichment(df_subset,df_dispersion,"es_bin","Effect size","AF_bin", "odds_ratio", 
                                        f"{dataset}, track {track_num}", 
                                        f"Plots_enrichment/match_{match}_af_enrichment_{dataset}_track{track_num}.pdf",
                                        color_mapping,
                                        match=match)
        df_subset.drop(columns=["es_bin"],inplace=True)








#----------
# Archived
#----------

# for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
#     df_subset=df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
#     df_subset["TF_type"]=get_contextual_tf_type(df_subset,dataset)
#     # plot histogram
#     sns.histplot(data=df_subset,x="log10_AF",hue="TF_type",bins=100)
#     plt.title(f"Allele frequency of SNVs overlapping REs of {dataset}")
#     plt.yscale("log")
#     plt.savefig(f"Plots/allele_frequency_distribution_{dataset}.pdf")
#     plt.close()
#     # plot percentage
#     sns.histplot(data=df_subset,x="log10_AF",hue="TF_type",bins=100,multiple="fill")
#     plt.title(f"Allele frequency of SNVs overlapping REs of {dataset}")
#     plt.savefig(f"Plots_maf/allele_frequency_distribution_fill_{dataset}.pdf")
#     plt.close()




# # for each dataset
# for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
#     df_subset=df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
#     df_subset["TF_type"]=get_contextual_tf_type(df_subset,dataset)
#     df_es_binned = df_subset.loc[:, ["AF_bin", "TF_type"]].copy().groupby(["AF_bin", "TF_type"]).size().unstack()
#     df_es_binned = calc_or(df_es_binned)
#     plot_or(transform_data_for_plotting(df_es_binned, "AF_bin"),
#             "AF_bin", "odds_ratio",dataset,{'sub': "#1f77b4", 'super': '#ff7f0e'},
#             out_name=f"Plots/af_enrichment_{dataset}.pdf",
#             ax=None)





# for cell_type in ["hepg2","k562"]:
#     df_subset=df_tfbs_maf_es[df_tfbs_maf_es["dataset"].str.contains(cell_type)].reset_index(drop=True)
#     if cell_type=="hepg2":
#         track_list=[0,2,4,6]
#     else:
#         track_list=[1,3,5,7]
#     for track_num in track_list:
#         df_subset=bin_and_label(df_subset,f"track_{track_num}", [-np.inf, -0.2, -0.1, 0, 0.1, 0.2, np.inf])
#         for i, af_bin in enumerate(df_es_binned['AF_bin'].unique()):
#             logger.info(f"Processing track {track_num}, allele frequency strata: {af_bin}")
#             df_es_binned=df_subset[df_subset["AF_bin"]==af_bin].reset_index(drop=True)
#             df_es_binned=df_es_binned.loc[:,["Bin","TF_type"]].copy().groupby(["Bin","TF_type"]).size().unstack()
#             df_es_binned=calc_or(df_es_binned)
#             df=transform_data_for_plotting(df_es_binned,"Bin") 
#             plot_or(df,"Bin","odds_ratio",
#                     f"Allele frequency strata: {af_bin}, track {track_num}",
#                     f"Plots_maf/es_enrichment_{af_bin}_track{track_num}.pdf",
#                     {'sub': "#1f77b4", 'super': '#ff7f0e'})



# def plot_contextual_sub_vs_super(df_orig, value_col):
#     df=df_orig.copy()
#     df[value_col]=df[value_col].abs()
#     for dataset in df["dataset"].unique():
#         df_subset = df[df["dataset"] == dataset].reset_index(drop=True).copy()
#         _, p_value = stats.mannwhitneyu(df_subset[df_subset["tf_property"] == "sub"][value_col],
#                                              df_subset[df_subset["tf_property"] == "super"][value_col])
#         sub_median = df_subset[df_subset["tf_property"] == "sub"][value_col].median()
#         super_median = df_subset[df_subset["tf_property"] == "super"][value_col].median()
#         larger_median_group = "sub" if sub_median > super_median else "super"
#         sns.kdeplot(data=df_subset, x=value_col, hue="tf_property", common_norm=False)
#         plt.title(f"{value_col} {dataset}")
#         plt.annotate(f"p={p_value:.1e}", xy=(0.5, 0.5), xycoords="axes fraction")
#         plt.annotate(f"larger median: {larger_median_group}", xy=(0.5, 0.4), xycoords="axes fraction")
#         plt.savefig(f"Plots_maf/contextual_sub_vs_super_{value_col}_{dataset}.pdf")
#         plt.close()


# df_tfbs_maf_es=get_additivity_profile(df_tfbs_maf_es)



# plot_contextual_sub_vs_super(df_tfbs_maf, "log10_AF")
# for i in range(8):
#     plot_contextual_sub_vs_super(df_tfbs_maf_es, f"track_{i}")


