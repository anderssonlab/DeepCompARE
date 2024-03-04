
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_and_label, replace_inf, get_minimum_positive, calc_odds_ratio




#--------------------------------------------------
# Functions
#--------------------------------------------------

def bin_variant(df,track_name):
    df=bin_and_label(df, "AF", [0,0.01,0.05,1])
    df["variant_type"]=["rare" if i=="0-0.01" else "low" if i=="0.01-0.05" else "common" for i in df["Bin"]]
    df.drop(columns=["Bin"],inplace=True)
    df=bin_and_label(df, track_name, [-np.inf,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,np.inf])
    df=df.loc[:,["variant_type","Bin"]].copy()
    df=df.groupby(["Bin","variant_type"]).size().unstack()
    return df



def plot_or(df_plot,title,out_name):
    df_plot["Predicted_effect_size"]=df_plot.index
    df_plot=df_plot.reset_index()
    df_plot=pd.melt(df_plot,id_vars=["Predicted_effect_size"],value_vars=["or_rare","or_low","or_common"],var_name="variant_type",value_name="odds_ratio")
    # in column variant_type, replace or_rare with rare, or_low with low, or_common with common
    df_plot["variant_type"]=df_plot["variant_type"].str.replace("or_","")
    # cat plot, use dashed line to connect the same variant type
    sns.catplot(data=df_plot,x="Predicted_effect_size",y="odds_ratio",hue="variant_type",kind="point",dodge=True,linestyles="--")    
    plt.title(title)
    plt.xlabel("Predicted effect size")
    plt.ylabel("Odds ratio")
    plt.xticks(rotation=45)
    plt.subplots_adjust(top=0.9,bottom=0.2)
    plt.savefig(out_name)
    plt.close()




# Odds ratio calculation 1: odds --> normalize column
def calc_or1(df):
    # calculate odds
    df["odds_rare"]=df["rare"]/(df["low"]+df["common"])
    df["odds_low"]=df["low"]/(df["rare"]+df["common"])
    df["odds_common"]=df["common"]/(df["rare"]+df["low"])
    # replace inf with second maximum unique value of the column * 1.5
    df["odds_low"]=replace_inf(df["odds_low"])
    df["odds_rare"]=replace_inf(df["odds_rare"])
    df["odds_common"]=replace_inf(df["odds_common"])
    # calculate odds ratio
    df["or_rare"]=df["odds_rare"]/get_minimum_positive(df["odds_rare"])
    df["or_low"]=df["odds_low"]/get_minimum_positive(df["odds_low"])
    df["or_common"]=df["odds_common"]/get_minimum_positive(df["odds_common"])
    return df

# Odds ratio calculation 2: canonical
def calc_or2(df):
    df_res=df.copy()
    df["or_rare"]=df_res.apply(lambda x: calc_odds_ratio(df_res,x.name,"rare"),axis=1)
    df["or_low"]=df_res.apply(lambda x: calc_odds_ratio(df_res,x.name,"low"),axis=1)
    df["or_common"]=df_res.apply(lambda x: calc_odds_ratio(df_res,x.name,"common"),axis=1)
    return df


#--------------------------------------------------
# Analysis 
#--------------------------------------------------

# Definition:
# Rare: AF<0.01
# Low: 0.01<=AF<0.05
# Common: AF>0.05

for file in ["enhancers_k562","enhancers_hepg2","promoters_k562","promoters_hepg2"]:
    for track_num in range(8):
        logger.info(f"Processing {file}, track {track_num}")
        df=pd.read_csv(f"maf_with_effect_size_{file}.csv",header=None)
        df.columns=["AF"]+["track_"+str(i) for i in range(16)]
        df=bin_variant(df,f"track_{track_num}")
        df=df[df.sum(axis=1)>10]
        df=calc_or2(df)
        plot_or(df,
                f"Odds ratio ({file}, track {track_num})",
                f"or_{file}_track{track_num}.png")
        logger.info(f"Done {file}, track {track_num}")

# nohup python3 plot_res.py > plot_res.out &

