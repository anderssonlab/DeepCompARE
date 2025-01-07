import pandas as pd
import numpy as np
from loguru import logger
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from stat_tests import bin_and_label, calc_or_by_various_thresholds

#-------------------
# Helper functions
#-------------------



# def bin_and_label(df, column_name, bin_edges):
#     bin_edges = sorted(set(bin_edges))
#     labels = [f"{bin_edges[i]} - {bin_edges[i+1]}" for i in range(len(bin_edges)-1)]
#     df[f"{column_name}_bin"] = pd.cut(df[column_name], bins=bin_edges, labels=labels, include_lowest=True)
#     df[f"{column_name}_bin"] = df[f"{column_name}_bin"].astype(pd.CategoricalDtype(categories=labels, ordered=True))
#     return df


# from scipy.stats import fisher_exact

# def odds_ratio_one_bin(df,bin_name):
#     """
#     Calculate odds ratio for a given bin
#     """
#     A=df.loc[bin_name,"True"]
#     B=df.loc[bin_name,"False"]
#     C=df.loc[:,"True"].sum()-A
#     D=df.loc[:,"False"].sum()-B
#     table = np.array([[A, B], [C, D]])
#     odds_ratio, p = fisher_exact(table)
#     se_log_or=np.sqrt(1/A+1/B+1/C+1/D)
#     ci_lower=np.exp(np.log(odds_ratio)-1.96*se_log_or)
#     ci_upper=np.exp(np.log(odds_ratio)+1.96*se_log_or)
#     return odds_ratio, p, ci_lower, ci_upper


# def odds_ratio_one_df(df,suffix):
#     """
#     Iterate over all bins (rows) and calculate odds ratio for each bin
#     """
#     if df.shape[0]==0:
#         logger.warning("Empty data frame")
#         return None,None
#     df=df.loc[(df["True"]>0) & (df["False"]>0)]
#     or_list = []
#     pval_list = []
#     ci_low_list = []
#     ci_high_list = []
#     for index, _ in df.iterrows():
#         or_value, pval_value, ci_low_value, ci_high_value = odds_ratio_one_bin(df, index)
#         or_list.append(or_value)
#         pval_list.append(pval_value)
#         ci_low_list.append(ci_low_value)
#         ci_high_list.append(ci_high_value)
#     df_res=pd.DataFrame({f"or":or_list,
#                          f"pval":pval_list,
#                          f"ci_low":ci_low_list,
#                          f"ci_high":ci_high_list},
#                         index=df.index)
#     df_res["threshold"]=suffix
#     return df_res



# def calc_or_by_various_thresholds(df,threshold_col,threshold_list,operation,bin_col):
#     df_res=pd.DataFrame()
#     n_list=[]
#     for threshold in threshold_list:
#         if operation=="larger":
#             df["target"]=(df[threshold_col]>threshold)
#         if operation=="smaller":
#             df["target"]=(df[threshold_col]<threshold)
#         n_list.append(df["target"].sum())
#         df["target"]=df["target"].astype(str)
#         df_plot=df.groupby([bin_col,"target"]).size().unstack()
#         df_temp=odds_ratio_one_df(df_plot,threshold)
#         df_res=pd.concat([df_res,df_temp])
#     return df_res,n_list








def plot_or(df,color_mapping,n_mapping,operation,title,out_path):
    """
    df: output of odds_ratio_one_df, contains 'or', 'pval', 'ci_low', 'ci_high', 'threshold'
    """
    plt.figure(figsize=(8, 6))
    # determine transparency
    df["alphas"] = df["pval"].apply(lambda x: 1.0 if x < 0.001 else (0.1 if x > 0.05 else 0.5))
    # x is the order of bin in categorical
    df["bin"]=df.index
    df["x"]=df["bin"].cat.codes
    jitter=0
    for threshold, df_subset in df.groupby(["threshold"]):
        plt.scatter(df_subset["x"]+jitter, 
                    df_subset["or"], 
                    color=color_mapping[threshold], 
                    alpha=df_subset['alphas'])
        for _, row in df_subset.iterrows():
            plt.errorbar(
                row["x"] + jitter, 
                row["or"], 
                yerr=[[row["or"] - row["ci_low"]], [row["ci_high"] - row["or"]]],
                fmt='o', 
                color=color_mapping[threshold],
                capsize=0, 
                alpha=row["alphas"],  # Set transparency for error bars
                markeredgewidth=1)
        if operation=="larger":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"allele frequency > {threshold} (n={n_mapping[threshold]})"
            )
        if operation=="smaller":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"allele frequency < {threshold} (n={n_mapping[threshold]})"
            )
        jitter+=0.1
    plt.xlabel("Predicted loss-of-function effect size")
    plt.ylabel("Odds ratio")
    plt.title(title)
    plt.axhline(y=1, color='black', linestyle=':')
    plt.legend()
    # change the ticks back to the original labels: bin edges
    plt.xticks(df["x"], df["bin"],rotation=45)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()





#----------------------
# Read data
#----------------------
df_promoters_hepg2=pd.read_csv("Pd1_maf_with_effect_size/maf_with_effect_size_promoters_hepg2.csv",header=None,index_col=0)
df_promoters_hepg2["dataset"]="promoters_hepg2"
df_promoters_k562=pd.read_csv("Pd1_maf_with_effect_size/maf_with_effect_size_promoters_k562.csv",header=None,index_col=0)
df_promoters_k562["dataset"]="promoters_k562"
df_enhancers_hepg2=pd.read_csv("Pd1_maf_with_effect_size/maf_with_effect_size_enhancers_hepg2.csv",header=None,index_col=0)
df_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_enhancers_k562=pd.read_csv("Pd1_maf_with_effect_size/maf_with_effect_size_enhancers_k562.csv",header=None,index_col=0)
df_enhancers_k562["dataset"]="enhancers_k562"



df_orig=pd.concat([df_promoters_hepg2,df_promoters_k562,df_enhancers_hepg2,df_enhancers_k562]).reset_index(drop=True)
df_orig.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+["track_"+str(i) for i in range(16)]+["dataset"]
df_orig["cell_type"]=df_orig["dataset"].str.split("_").str[1]


threshold_list = [0.0001,0.001,0.01,0.05]


dataset="promoters_hepg2"
track_num=0

for dataset in ["promoters_hepg2","promoters_k562","enhancers_hepg2","enhancers_k562"]:
    for track_num in get_track_num(dataset):
        logger.info(f"Processing {dataset}, track {track_num}")
        df=df_orig[df_orig["dataset"]==dataset].reset_index(drop=True)
        df=df[df[f"track_{track_num}"]<0].reset_index(drop=True)
        df[f"track_{track_num}"]=df[f"track_{track_num}"].abs()
        df=bin_and_label(df,f"track_{track_num}",[0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])
        df_or,n_list=calc_or_by_various_thresholds(df,"AF",threshold_list,"larger",f"track_{track_num}_bin")
        colors_list = ['#808080', '#c9d800', '#27ae60', 'blue']
        color_mapping = dict(zip(threshold_list, colors_list))
        n_mapping = dict(zip(threshold_list, n_list))
        plot_or(df_or,color_mapping,n_mapping,"larger",f"Odds ratio ({dataset}, track {track_num})",f"or_{dataset}_track{track_num}.pdf")


# nohup python3 plot_or.py > plot_or.out &
