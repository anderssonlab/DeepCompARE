import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from stat_tests import bin_and_label, calc_or_by_various_thresholds


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
#-------------------
# Helper functions
#-------------------
def read_file(file_suffix):
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv',index_col=0)
    df["max_af"] = df["gnomad_af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["max_241way"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    return df



def plot_or_af(df,color_mapping,n_mapping,operation,title,out_path):
    """
    df: output of odds_ratio_one_df, contains 'or', 'pval', 'ci_low', 'ci_high', 'threshold'
    """
    plt.figure(figsize=(2.65,2.15))
    # thinner frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
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
                    alpha=df_subset['alphas'],
                    s=7)
        for _, row in df_subset.iterrows():
            plt.errorbar(
                row["x"] + jitter, 
                row["or"], 
                yerr=[[row["or"] - row["ci_low"]], [row["ci_high"] - row["or"]]],
                fmt='none', 
                color=color_mapping[threshold],
                capsize=0, 
                alpha=row["alphas"],  # Set transparency for error bars
                markeredgewidth=0.5,
                elinewidth=1)
        if operation=="larger":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"max(AF)>{threshold} (n={n_mapping[threshold]})",
                s=10)
        if operation=="smaller":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"max(AF)<{threshold} (n={n_mapping[threshold]})",
                s=10)
        jitter+=0.15
    plt.xlabel("Motif ISA score",fontsize=7)
    plt.ylabel("Odds ratio",fontsize=7)
    plt.axhline(y=1, color='black', linestyle=':',linewidth=0.5)
    plt.legend(fontsize=5)
    # change the ticks back to the original labels: bin edges
    plt.xticks(df["x"], df["bin"],rotation=30,fontsize=5)
    plt.yticks(fontsize=5)
    # small margin
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(out_path,dpi=300)
    plt.close()










def plot_or_phylop(df,color_mapping,n_mapping,operation,title,out_path):
    """
    df: output of odds_ratio_one_df, contains 'or', 'pval', 'ci_low', 'ci_high', 'threshold'
    """
    plt.figure(figsize=(2.65, 2.15))
    # thinner frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
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
                    alpha=df_subset['alphas'],
                    s=7)
        for _, row in df_subset.iterrows():
            plt.errorbar(
                row["x"] + jitter, 
                row["or"], 
                yerr=[[row["or"] - row["ci_low"]], [row["ci_high"] - row["or"]]],
                fmt='none',
                color=color_mapping[threshold],
                capsize=0, 
                alpha=row["alphas"],  # Set transparency for error bars
                markeredgewidth=0.5,
                elinewidth=1)
        if operation=="larger":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"max(phylop)>{threshold} (n={n_mapping[threshold]})",s=10
            )
        if operation=="smaller":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"max(phylop)<{threshold} (n={n_mapping[threshold]})",s=10
            )
        jitter+=0.15
    plt.xlabel("Motif ISA score",fontsize=7)
    plt.ylabel("Odds ratio",fontsize=7)
    plt.axhline(y=1, color='black', linestyle=':',linewidth=0.5)
    # legend at top left
    plt.legend(fontsize=5,loc='upper left')
    # change the ticks back to the original labels: bin edges
    plt.xticks(df["x"], df["bin"],rotation=30,fontsize=5)
    plt.yticks(fontsize=5)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()





#--------------------------------------------------------------
# Analysis 1: effect size v.s. allele frequency
# Goal: reproduce SNP level enrichment
#--------------------------------------------------------------


# Read in data
df_promoter_hepg2=read_file("promoters_hepg2")
df_enhancer_hepg2=read_file("enhancers_hepg2")
df_promoter_k562=read_file("promoters_k562")
df_enhancer_k562=read_file("enhancers_k562")

df_orig=pd.concat([df_promoter_hepg2,df_enhancer_hepg2,df_promoter_k562,df_enhancer_k562], axis=0)
df_orig["cell_type"]=df_orig["dataset"].apply(lambda x: x.split("_")[1])




colors_list = ['#808080', '#B5BD00', '#27ae60', '#1E90FF']
# plot max_af
threshold_list = [0.0001,0.001,0.01,0.05]
for dataset in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
    for track_num in get_track_num(dataset):
        logger.info(f"Processing {dataset}, track {track_num}")
        df=df_orig[df_orig["dataset"]==dataset].reset_index(drop=True)
        df=bin_and_label(df,f"isa_track{track_num}",[-np.inf, -0.8, -0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])
        # for common variants
        df_or,n_list=calc_or_by_various_thresholds(df,"max_af",threshold_list,"larger",f"isa_track{track_num}_bin")
        color_mapping = dict(zip(threshold_list, colors_list))
        n_mapping = dict(zip(threshold_list, n_list))
        plot_or_af(df_or,color_mapping,n_mapping,"larger",f"Odds ratio ({dataset}, track {track_num})",f"Plots_af/af_or_{dataset}_track{track_num}.pdf")



# colors_list = ["#A9A9A9", "#FFC0CB", "#C7B6EA", "#9400D3"] pink to purple
colors_list = ["#808080", '#FFD700', '#FFA500', '#D98880']
# plot phylop
threshold_list = [0,1,2,3]
for dataset in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
    for track_num in get_track_num(dataset):
        logger.info(f"Processing {dataset}, track {track_num}")
        df=df_orig[df_orig["dataset"]==dataset].reset_index(drop=True)
        df=bin_and_label(df,f"isa_track{track_num}",[-np.inf, -0.8, -0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])
        df_or,n_list=calc_or_by_various_thresholds(df,"max_241way", threshold_list,"larger",f"isa_track{track_num}_bin")
        color_mapping = dict(zip(threshold_list, colors_list))
        n_mapping = dict(zip(threshold_list, n_list))
        plot_or_phylop(df_or,color_mapping,n_mapping,"larger",f"Odds ratio ({dataset}, track {track_num})",f"Plots_phylop/phylop_or_{dataset}_track{track_num}.pdf")



# nohup python3 plot_or.py > plot_or.out &