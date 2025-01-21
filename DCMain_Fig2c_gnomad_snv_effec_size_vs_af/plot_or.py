import pandas as pd
import numpy as np
from loguru import logger
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from stat_tests import bin_and_label, calc_or_by_various_thresholds


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42

def plot_or(df,color_mapping,n_mapping,operation,title,out_path):
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
                label=f"AF > {threshold} (n={n_mapping[threshold]})",s=10)
        if operation=="smaller":
            plt.scatter(
                [], [],  # Invisible data
                color=color_mapping[threshold], 
                label=f"AF < {threshold} (n={n_mapping[threshold]})",s=10)
        jitter+=0.15
    plt.xlabel("Predicted loss-of-function effect size", fontsize=7)
    plt.ylabel("Odds ratio", fontsize=7)
    # plt.title(title, fontsize=7)
    plt.axhline(y=1, color='black', linestyle=':', linewidth=0.5)
    plt.legend(fontsize=5,loc='lower left')
    # change the ticks back to the original labels: bin edges
    plt.xticks(df["x"], df["bin"],rotation=30,fontsize=5)
    plt.yticks(fontsize=5)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
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




for dataset in ["promoters_hepg2","promoters_k562","enhancers_hepg2","enhancers_k562"]:
    for track_num in get_track_num(dataset):
        logger.info(f"Processing {dataset}, track {track_num}")
        df=df_orig[df_orig["dataset"]==dataset].reset_index(drop=True)
        df=df[df[f"track_{track_num}"]<0].reset_index(drop=True)
        df[f"track_{track_num}"]=df[f"track_{track_num}"].abs()
        df=bin_and_label(df,f"track_{track_num}",[0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])
        df_or,n_list=calc_or_by_various_thresholds(df,"AF",threshold_list,"larger",f"track_{track_num}_bin")
        colors_list = ['#808080', '#B5BD00', '#27ae60', '#1E90FF']
        color_mapping = dict(zip(threshold_list, colors_list))
        n_mapping = dict(zip(threshold_list, n_list))
        plot_or(df_or,color_mapping,n_mapping,"larger",f"Odds ratio ({dataset}, track {track_num})",f"or_{dataset}_track{track_num}.pdf")


# nohup python3 plot_or.py > plot_or.out &
