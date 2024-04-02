
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_and_label, calc_odds_ratio




#--------------------------------------------------
# Functions
#--------------------------------------------------

def bin_variant(df,track_name):
    df=bin_and_label(df, "AF", [0,0.01,0.05,1])
    df["variant_type"]=["rare" if i=="0-0.01" else "low" if i=="0.01-0.05" else "common" for i in df["Bin"]]
    df.drop(columns=["Bin"],inplace=True)
    df=bin_and_label(df, track_name, [-np.inf,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,np.inf])
    df=df.loc[:,["variant_type","Bin"]].copy()
    df=df.groupby(["Bin","variant_type"]).size().unstack()
    return df


def calc_or(df):
    df_res=df.copy()
    df[['or_rare', 'pval_rare']]=df_res.apply(lambda x: calc_odds_ratio(df_res,x.name,"rare"),axis=1).to_list()
    df[['or_low', 'pval_low']]=df_res.apply(lambda x: calc_odds_ratio(df_res,x.name,"low"),axis=1).to_list()
    df[['or_common', 'pval_common']]=df_res.apply(lambda x: calc_odds_ratio(df_res,x.name,"common"),axis=1).to_list()
    return df



def plot_or(df_plot,title,out_name):
    df_plot["Predicted_effect_size"]=df_plot.index
    df_plot=df_plot.reset_index()
    df_or=pd.melt(df_plot,id_vars=["Predicted_effect_size"],value_vars=["or_rare","or_low","or_common"],var_name="variant_type",value_name="odds_ratio")
    df_pval=pd.melt(df_plot,id_vars=["Predicted_effect_size"],value_vars=["pval_rare","pval_low","pval_common"],var_name="variant_type",value_name="pval")
    df_or["variant_type"]=df_or["variant_type"].str.replace("or_","")
    df_pval["variant_type"]=df_pval["variant_type"].str.replace("pval_","")
    df_plot=pd.merge(df_or,df_pval,on=["Predicted_effect_size","variant_type"])
    df_plot["alphas"]=df_plot["pval"].apply(lambda x: 1.0 if x < 0.001 else (0.1 if x > 0.05 else 0.5))

    plt.figure(figsize=(4, 4))
    color_mapping = {'rare': "#1f77b4", 'low': '#ff7f0e', 'common': '#2ca02c'}

    # Plot each point individually
    for variant, df_subset in df_plot.groupby('variant_type'):
        # TODO: choose whether to plot low
        if variant == 'low':
            continue
    # Sort the subset for consistent plotting
        # Plotting the lines
        plt.plot(df_subset['Predicted_effect_size'], df_subset['odds_ratio'], '--', color=color_mapping[variant], label=variant)
        # Plotting the dots
        plt.scatter(df_subset['Predicted_effect_size'], df_subset['odds_ratio'], color=color_mapping[variant], alpha=df_subset['alphas'])
        
    plt.title(title)
    plt.xlabel("Predicted effect size")
    plt.ylabel("Odds ratio")
    plt.xticks(rotation=45)
    plt.legend(title='Variant Type')
    plt.tight_layout()
    plt.savefig(out_name)
    plt.close()




#--------------------------------------------------
# Analysis 
#--------------------------------------------------

# Definition:
# Rare: AF<0.01
# Low: 0.01<=AF<0.05
# Common: AF>0.05


def get_track_num(file):
    if "k562" in file:
        return list(range(1,16,2))
    elif "hepg2" in file:
        return list(range(0,16,2))
    return None


for file in ["enhancers_k562","enhancers_hepg2","promoters_k562","promoters_hepg2"]:
    for track_num in get_track_num(file):
        logger.info(f"Processing {file}, track {track_num}")
        df=pd.read_csv(f"maf_with_effect_size_{file}.csv",header=None)
        df.columns=["AF"]+["track_"+str(i) for i in range(16)]
        df=bin_variant(df,f"track_{track_num}")
        df=df[df.sum(axis=1)>10]
        df=calc_or(df)
        df.to_csv("Or_pval_tables/"+f"{file}_track{track_num}.csv",index=False)
        # TODO: change directory whether to plot low
        plot_or(df,
                f"Odds ratio ({file}, track {track_num})",
                f"Plots_wo_low/or_{file}_track{track_num}.pdf")
        logger.info(f"Done {file}, track {track_num}")

# nohup python3 plot_res.py > plot_res.out &
