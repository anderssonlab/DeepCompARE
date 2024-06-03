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
    df=bin_and_label(df, track_name, [-np.inf,-0.2,0.2,np.inf])
    df["variant_type"]=["large_negative" if i=="-inf - -0.2" else "small" if i=="-0.2 - 0.2" else "large_positive" for i in df["Bin"]]
    df.drop(columns=["Bin"],inplace=True)
    df=bin_and_label(df, "log10_AF", [-np.inf,-4,-2,0])
    return df.loc[:,["variant_type","Bin"]].copy().groupby(["Bin","variant_type"]).size().unstack()

def calc_or(df):
    df_res = df.copy()
    categories = df.columns.tolist()
    for category in categories:
        df[[f'or_{category}', f'pval_{category}', f'ci_low_{category}', f'ci_high_{category}']] = df_res.apply(lambda x: calc_odds_ratio(df_res, x.name, category), axis=1).to_list()
    return df

def transform_data_for_plotting(df_plot):
    df_plot["log10(Allele frequency)"]=df_plot.index
    df_plot=df_plot.reset_index(drop=True)
    # melt columns
    or_columns=[i for i in df_plot.columns if i.startswith("or")]
    pval_columns=[i for i in df_plot.columns if i.startswith("pval")]
    ci_low_columns=[i for i in df_plot.columns if i.startswith("ci_low")]
    ci_high_columns=[i for i in df_plot.columns if i.startswith("ci_high")]
    df_or=pd.melt(df_plot,id_vars=["log10(Allele frequency)"],value_vars=or_columns,var_name="variant_type",value_name="odds_ratio")
    df_pval=pd.melt(df_plot,id_vars=["log10(Allele frequency)"],value_vars=pval_columns,var_name="variant_type",value_name="pval")
    df_ci_low=pd.melt(df_plot,id_vars=["log10(Allele frequency)"],value_vars=ci_low_columns,var_name="variant_type",value_name="ci_low")
    df_ci_high=pd.melt(df_plot,id_vars=["log10(Allele frequency)"],value_vars=ci_high_columns,var_name="variant_type",value_name="ci_high")
    # rename columns
    df_or["variant_type"]=df_or["variant_type"].str.replace("or_","")
    df_pval["variant_type"]=df_pval["variant_type"].str.replace("pval_","")
    df_ci_low["variant_type"]=df_ci_low["variant_type"].str.replace("ci_low_","")
    df_ci_high["variant_type"]=df_ci_high["variant_type"].str.replace("ci_high_","")
    # merge data frames
    df_plot=pd.merge(df_or,df_pval,on=["log10(Allele frequency)","variant_type"])
    df_plot=pd.merge(df_plot,df_ci_low,on=["log10(Allele frequency)","variant_type"])
    df_plot=pd.merge(df_plot,df_ci_high,on=["log10(Allele frequency)","variant_type"])
    # determine transparency
    df_plot["alphas"]=df_plot["pval"].apply(lambda x: 1.0 if x < 0.001 else (0.1 if x > 0.05 else 0.5))
    return df_plot

def convert_interval_to_midpoint(interval_str):
    start, end = map(float, interval_str.split(' - '))
    midpoint = (start + end) / 2
    return midpoint

def plot_or(df_plot,title,out_name):
    df_plot = transform_data_for_plotting(df_plot)
    plt.figure(figsize=(5, 5)) 
    color_mapping = {
        "large_negative": "#1f77b4",
        "medium_negative": "#ff7f0e",
        "small": "#2ca02c",
        "medium_positive": "#d62728",
        "large_positive": "#9467bd"
    }
    jitter_strength = 0.08 
    original_x_labels = {}
    for variant, df_subset in df_plot.groupby('variant_type'):
        x_values = df_subset['log10(Allele frequency)']
        y_values = df_subset['odds_ratio']
        x_values_numeric = np.float_(x_values.apply(convert_interval_to_midpoint))
        x_values_numeric[0] = x_values_numeric[1]-2
        for original_label, midpoint in zip(x_values, x_values_numeric):
            original_x_labels[midpoint] = original_label
        jittered_x_values = x_values_numeric + jitter_strength 
        jitter_strength+=0.08
        
        plt.plot(jittered_x_values, y_values, '--', color=color_mapping[variant], label=variant)
        plt.scatter(jittered_x_values, y_values, color=color_mapping[variant], alpha=df_subset['alphas'])
        plt.errorbar(jittered_x_values, 
                    y_values, 
                    yerr=[y_values - df_subset['ci_low'], df_subset['ci_high'] - y_values],
                    fmt='none', 
                    color=color_mapping[variant],
                    capsize=3, markeredgewidth=1)
    plt.title(title)
    plt.axhline(y=1, color='black', linestyle=':')
    plt.xlabel("log10(AF)")
    plt.ylabel("Odds ratio")
    midpoints = sorted(original_x_labels.keys())
    original_labels = [original_x_labels[midpoint] for midpoint in midpoints]
    plt.xticks(midpoints, original_labels, rotation=45)
    plt.legend(title='Variant Type')
    plt.tight_layout()
    plt.savefig(out_name)
    plt.close()


#--------------------------------------------------
# Analysis 
#--------------------------------------------------

def get_track_num(file):
    if "k562" in file:
        return list(range(1,8,2))
    elif "hepg2" in file:
        return list(range(0,8,2))
    return None


for file in ["enhancers_k562","enhancers_hepg2","promoters_k562","promoters_hepg2"]:
    for track_num in get_track_num(file):
        logger.info(f"Processing {file}, track {track_num}")
        df=pd.read_csv(f"maf_with_effect_size_{file}.csv",header=None,index_col=0)
        df.reset_index(drop=True,inplace=True)
        df.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
        df["log10_AF"]=np.log10(df["AF"])
        df=bin_variant(df,f"track_{track_num}")
        df=df[df.sum(axis=1)>10]
        df=df.loc[:, (df != 0).all(axis=0)]
        df=calc_or(df)
        df.to_csv("Or_pval_tables/"+f"{file}_track{track_num}.csv",index=False)
        plot_or(df,
                f"Odds ratio ({file}, track {track_num})",
                f"Plot_or_reverse/or_{file}_track{track_num}.pdf")
        logger.info(f"Done {file}, track {track_num}")

# nohup python3 plot_or_reverse.py > plot_or_reverse.out &



#------------
# Debugging
#------------

# file="enhancers_k562"
# track_num=3
# df=pd.read_csv(f"maf_with_effect_size_{file}.csv",header=None,index_col=0)
# df.reset_index(drop=True,inplace=True)
# df.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
# df["log10_AF"]=np.log10(df["AF"])
# df=bin_variant(df,f"track_{track_num}")
# df=df[df.sum(axis=1)>10]
# # remove column if any value in the column is zero
# df=df.loc[:, (df != 0).all(axis=0)]
# df=calc_or(df)
# plot_or(df,f"Odds ratio ({file}, track {track_num})",f"Plot_or_reverse/or_{file}_track{track_num}.pdf")