import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from loguru import logger



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import calc_or, plot_or
from utils import split_dimer




#-----------------------------------------------------
# helper dataset
#-----------------------------------------------------

# for match_cell_type_specificity()
df_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Cooperativity_effects_tf_level/TFs.dispersionEstimates.hepG2.tab",sep="\t")
df_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Cooperativity_effects_tf_level/TFs.dispersionEstimates.k562.tab",sep="\t")
df_dispersion=pd.concat([df_dispersion_hepg2,df_dispersion_k562],axis=0).reset_index(drop=True)
df_dispersion=df_dispersion.drop_duplicates().reset_index(drop=True)
df_dispersion["gini_rank"]=df_dispersion["gini"].rank(ascending=False)

#-----------------------------------------------------
# helper functions 
#-----------------------------------------------------

color_mapping = {
    'redundant': "#1f77b4", 
    'codependent': '#ff7f0e'
    }


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


def match_cell_type_specificity(df,df_dispersion,num_closest=1):
    redundant_tfs = df[df["TF_type"] == "redundant"]["motif"].unique().tolist()
    codependent_tfs = df[df["TF_type"] == "codependent"]["motif"].unique().tolist()
    redundant_ranks=df_dispersion[df_dispersion["gene"].isin(redundant_tfs)]["gini_rank"].to_list()
    codependent_ranks=df_dispersion[df_dispersion["gene"].isin(codependent_tfs)]["gini_rank"].to_list()
    matched_ranks=find_closest_numbers(redundant_ranks,codependent_ranks,num_closest)
    matched_tfs=df_dispersion[df_dispersion["gini_rank"].isin(matched_ranks)]["gene"].to_list()
    df_subset=df[df["motif"].isin(redundant_tfs+matched_tfs)].reset_index(drop=True).copy()
    return df_subset


def read_file(file_suffix):
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv')
    cols_remove=[col for col in df.columns if col.startswith('pred_orig')]
    df.drop(cols_remove, axis=1, inplace=True)
    df["max_af"] = df["af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["min_af"] = df["af"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["num_variants"] = df["af"].apply(lambda x: len(str(x).split(":")))
    df["num_common_variants"] = df["af"].apply(lambda x: sum([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["have_common_variant"]= (df["max_af"]>0.001)
    df["241way_max"] = df["241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["447way_max"] = df["447way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    return df




def add_tf_codependency(df):
    tfs_codependent=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_codependent.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_redundant.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df





df=read_file("promoters_hepg2")
df=add_tf_codependency(df)



def moving_thresh_odds_ratio(value_list,target_bool,thresh_list):
    """
    Calculate the odds ratio of the target_list when the value_list is greater than or equal to the threshold
    Args:
    value_list: list of values to be thresholded
    target_list: boolean list of target
    thresh_list: list of thresholds. The larger the threshold, the more likely the value is to be True
    """
    odds_ratios=[]
    for thresh in thresh_list:
        value_bool=[value>=thresh for value in value_list]
        df=pd.DataFrame({'value_bool': value_bool,'target_bool': target_bool})
        counts_TT=df[(df["target_bool"]==True) & (df["value_bool"]==True)].shape[0]
        counts_TF=df[(df["target_bool"]==True) & (df["value_bool"]==False)].shape[0]
        counts_FT=df[(df["target_bool"]==False) & (df["value_bool"]==True)].shape[0]
        counts_FF=df[(df["target_bool"]==False) & (df["value_bool"]==False)].shape[0]
        odds_true=(counts_TT+1)/(counts_TF+1)
        odds_false=(counts_FT+1)/(counts_FF+1)
        odds_ratios.append(odds_true/odds_false)
    return odds_ratios


def get_percentile_thresholds(df,value_col):
    """
    Get linearly spaced threshold
    target_tf_type: ensure that non-target TFs are included throughout the threshold range
    """
    df_nonrare=df[df["AF_bin"]!="rare"].reset_index(drop=True)
    percentiles=np.percentile(df_nonrare[value_col],list(range(5,98)))
    return percentiles



def moving_or(df,target_col,target_type,value_col,thresh_list=None):
    """
    target_col: AF_bin
    target_type: rare/common
    """
    if thresh_list is None:
        thresh_list=get_percentile_thresholds(df,value_col)
    df_redundant=df[df["TF_type"]=="redundant"].reset_index(drop=True)
    df_codependent=df[df["TF_type"]=="codependent"].reset_index(drop=True)
    or_redundant=moving_thresh_odds_ratio(df_redundant[value_col],df_redundant[target_col]==target_type,thresh_list)
    or_codependent=moving_thresh_odds_ratio(df_codependent[value_col],df_codependent[target_col]==target_type,thresh_list)
    df_plot_redundant=pd.DataFrame({"effect_size":thresh_list,"odds_ratio":or_redundant,"TF_type":"redundant"})
    df_plot_codependent=pd.DataFrame({"effect_size":thresh_list,"odds_ratio":or_codependent,"TF_type":"codependent"})
    df_plot=pd.concat([df_plot_redundant,df_plot_codependent],axis=0).reset_index(drop=True)
    return df_plot



def plot_rare_enrichment(df,cell_type,track_num):
    df = df.loc[:, ["TF_type","AF_bin",f"track_{track_num}"]].copy()
    df_plot=moving_or(df,"AF_bin","rare",f"track_{track_num}")
    sns.lineplot(data=df_plot,x="effect_size",y="odds_ratio",hue="TF_type") # ,style="value_type"
    plt.title(f"Enrichment of extremely rare variant, {cell_type}, track {track_num}")
    plt.savefig(f"Plots_rare_enrichment/rare_enrichment_{cell_type}_track{track_num}.pdf")
    plt.close()



for cell_type in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562","hepg2","k562"]:
    if "hepg2" in cell_type:
        track_list = [0, 2, 4, 6]
    else:
        track_list = [1, 3, 5, 7]
    df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"].str.contains(cell_type)].reset_index(drop=True)
    for track_num in track_list:
        logger.info(f"Processing cell type {cell_type}, track {track_num}")
        # subset to contain only negative effect size
        df_subset_neg=df_subset[df_subset[f"track_{track_num}"]<0].reset_index(drop=True)
        df_subset_neg[f"track_{track_num}"]=df_subset_neg[f"track_{track_num}"].abs()
        plot_rare_enrichment(df_subset_neg,cell_type,track_num)







# nohup python3 es_af_rare_enrichment.py > es_af_rare_enrichment.out &































#----------
# Archived
#----------


#--------------------------------------
# variant enrichment per bin
#--------------------------------------

# # TODO: summarize the enrichment sub in large effect size, in rare, low common variants

# for cell_type in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562","hepg2","k562"]:
#     if "hepg2" in cell_type:
#         track_list = [0, 2, 4, 6]
#     else:
#         track_list = [1, 3, 5, 7]
#     variant_type="rare"
#     df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"].str.contains(cell_type)].reset_index(drop=True)
#     df_subset = df_subset[df_subset["AF_bin"]==variant_type].reset_index(drop=True)
#     for track_num in track_list:
#         logger.info(f"Processing cell type {cell_type}, track {track_num}")
#         df_subset_neg=df_subset[df_subset[f"track_{track_num}"]<0].reset_index(drop=True)
#         df_subset_neg[f"track_{track_num}"]=df_subset_neg[f"track_{track_num}"].abs()
#         sns.kdeplot(data=df_subset_neg,x=f"track_{track_num}",hue="TF_type",common_norm=False,cumulative=True)
#         plt.title(f"cumulative distribution of effect size for {variant_type} variant\n{cell_type} track {track_num}")
#         plt.axhline(y=1, color="black",linestyle=':')
#         # add two vertical lines to show the 95 percentile effect size for sub and super TFs
#         plt.axvline(x=df_subset_neg[df_subset_neg["TF_type"]=="redundant"][f"track_{track_num}"].quantile(0.99),color=color_mapping["redundant"],linestyle='--')
#         plt.axvline(x=df_subset_neg[df_subset_neg["TF_type"]=="codependent"][f"track_{track_num}"].quantile(0.99),color=color_mapping["codependent"],linestyle='--')
#         plt.savefig(f"Plots_{variant_type}_distribution/cumulative_distribution_{variant_type}_{cell_type}_track_{track_num}.pdf")
#         plt.close()
        







# def get_or_for_plot(df,x_colname_to_bin,hue_colname="TF_type"):
#     bins=np.linspace(0,df[x_colname_to_bin].max(),10)
#     df= bin_and_label(df, f"track_{track_num}",bins,bin_names="midpoint",new_column_name="es_bin")
#     df = df.loc[:, ["es_bin", hue_colname]].copy().groupby(["es_bin", hue_colname]).size().unstack()
#     df = df.fillna(0)
#     df+=1
#     print(df)
#     df_or = calc_or(df, "es_bin", hue_colname, out_group="other")
#     return df_or


# def plot_stratified_enrichment(df,strata_colname,strata_prefix, x_colname_to_bin, sup_title, out_name, color_mapping):
#     strata_bins = df[strata_colname].unique()
#     fig, axs = plt.subplots(len(strata_bins), 1, figsize=(6, 8))
#     for i, strata_bin in enumerate(strata_bins):
#         logger.info(f"Processing {strata_prefix}: {strata_bin}")
#         df_strata_subset = df[df[strata_colname] == strata_bin].reset_index(drop=True)
#         if df_strata_subset.shape[0] < 10:
#             continue
#         # match super to sub
#         df_plot=get_or_for_plot(df_strata_subset, x_colname_to_bin)
#         plot_or(df_plot, "es_bin", "odds_ratio", "TF_type",
#                 f"{strata_prefix}: {strata_bin}", color_mapping, ax=axs[i])
#     axs[-1].set_xlabel(x_colname_to_bin)
#     fig.text(0.04, 0.5, "odds_ratio", va='center', rotation='vertical')
#     plt.suptitle(sup_title) 
#     plt.tight_layout()  
#     plt.xticks(rotation=45)
#     plt.subplots_adjust(bottom=0.15)
#     plt.subplots_adjust(left=0.15)
#     plt.savefig(out_name)
#     plt.close()
    


# # dataset="enhancers_hepg2"
# # track_num=0
# for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
#     df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
#     sub_tfs,super_tfs=get_sub_super_tfs(dataset)
#     matched_super_tfs=match_super_with_sub(sub_tfs,super_tfs,df_dispersion) # no matching, no expected results.
#     df_subset=df_subset[df_subset["motif"].isin(sub_tfs+matched_super_tfs)].reset_index(drop=True)
#     df_subset["TF_type"]=np.where(df_subset["motif"].isin(sub_tfs),"sub","super")
#     if "hepg2" in dataset:
#         track_list = [0, 2, 4, 6]
#         df_dispersion=df_dispersion_hepg2
#     else:
#         track_list = [1, 3, 5, 7]
#         df_dispersion=df_dispersion_k562
#     for track_num in track_list:
#         logger.info(f"Processing dataset {dataset}, track {track_num}")
#         df_subset_neg=df_subset[df_subset[f"track_{track_num}"]<0].reset_index(drop=True)
#         df_subset_neg[f"track_{track_num}"]=df_subset_neg[f"track_{track_num}"].abs()
#         logger.info(f"Processing dataset {dataset}, track {track_num}")
#         # stratify by AF_bin, plot es_enrichment
#         plot_stratified_enrichment(df_subset_neg,"AF_bin","Variant type",f"track_{track_num}", 
#                                     f"{dataset}, track {track_num}", 
#                                     f"Plots_enrichment/es_enrichment_{dataset}_track{track_num}.pdf",
#                                     color_mapping)
#         plot_stratified_enrichment(df_subset_neg,df_dispersion,"es_bin","Effect size","AF_bin", "odds_ratio", 
#                         f"{dataset}, track {track_num}", 
#                         f"Plots_enrichment/match_{match}_af_enrichment_{dataset}_track{track_num}.pdf",
#                         color_mapping,
#                         match=match)



# to compare between positives and negatives
# def moving_or_pos_neg(df,target_col,target_type,value_col):
#     df_neg=df[df[value_col]<0].reset_index(drop=True)
#     df_pos=df[df[value_col]>=0].reset_index(drop=True)
#     df_neg[value_col]=df_neg[value_col].abs()
#     df=pd.concat([df_neg,df_pos],axis=0).reset_index(drop=True)
#     thresh_list=get_percentile_thresholds(df,target_col,target_type,value_col)
#     df_plot_pos=moving_or_sub_super(df_pos,target_col,target_type,value_col,thresh_list)
#     df_plot_pos["value_type"]="positive"
#     df_plot_neg=moving_or_sub_super(df_neg,target_col,target_type,value_col,thresh_list)
#     df_plot_neg["value_type"]="negative"
#     df_plot=pd.concat([df_plot_pos,df_plot_neg],axis=0).reset_index(drop=True)
#     return df_plot


# def plot_rare_enrichment_pos_neg(df_orig,cell_type,track_num):
#     df = df_orig.loc[:, ["TF_type","AF_bin",f"track_{track_num}"]].copy()
#     df_plot=moving_or_pos_neg(df,"AF_bin","rare",f"track_{track_num}")
#     sns.lineplot(data=df_plot,x="effect_size",y="odds_ratio",hue="TF_type" ,style="value_type")
#     plt.title(f"Rare variant enrichment, {cell_type}, track {track_num}")
#     plt.savefig(f"Plots_rare_enrichment_ep/pos_neg_rare_enrichment_{cell_type}_track{track_num}.pdf")
#     plt.close()

# for cell_type in ["hepg2","k562"]:
#     if cell_type=="hepg2":
#         track_list = [0, 2, 4, 6]
#     else:
#         track_list = [1, 3, 5, 7]
#     for track_num in track_list:
#         logger.info(f"Processing cell type {cell_type}, track {track_num}")
#         df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"].str.contains(cell_type)].reset_index(drop=True)
#         df_subset["TF_type"]=get_contextual_tf_type(df_subset)
#         plot_rare_enrichment_pos_neg(df_subset,cell_type,track_num)













# method 1 percentile thresholds
# thresh_list=np.percentile(df_subset_neg_nonrare[f"track_{track_num}"],list(range(5,100)))
# df_map=pd.DataFrame({"percentile":list(range(5,100)),"value":thresh_list}).set_index("percentile")
# df_map.loc[100,"value"]=max(df_subset_neg_nonrare[f"track_{track_num}"])

# fig, ax1 = plt.subplots()
# sns.lineplot(data=df_plot,x="percentile",y="odds_ratio",hue="TF_type",ax=ax1)
# ax1.set_xticks([20, 40, 60, 80, 100])
# ax2 = ax1.twiny()
# # given percentile in primary_ticks, find the original value
# secondary_x_ticks = df_map.loc[[20, 40, 60, 80, 100],"value"]
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(ax1.get_xticks())
# # set_xticklabels to secondary_x_ticks, but in 1e-2 format
# ax2.set_xticklabels([f"{x:.2e}" for x in secondary_x_ticks])
# ax2.set_xlabel('Effect size')
# plt.savefig("odds_ratio_vs_effect_size.pdf")
# plt.close()





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


