import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_and_label



def get_track_num(file_suffix):
    if "hepg2" in file_suffix:
        return {"cage":0,"dhs":2,"starr":4,"sure":6}
    if "k562" in file_suffix:
        return {"cage":1,"dhs":3,"starr":5,"sure":7}


def read_file(file_suffix):
    df=pd.read_csv(f"maf_with_effect_size_{file_suffix}.csv",header=None,index_col=0)
    df.reset_index(drop=True,inplace=True)
    df.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
    df["log10_AF"]=np.log10(df["AF"])
    return df

def bin_file(df,track_num):
    df=bin_and_label(df, f"track_{track_num}", [-np.inf,-0.2,0.2,np.inf])
    df["variant_type"]=["large_negative" if i=="-inf - -0.2" else "small" if i=="-0.2 - 0.2" else "large_positive" for i in df["Bin"]]
    df.drop(columns=["Bin"],inplace=True)
    return df


def plot_file(file_suffix):
    df=read_file(file_suffix)
    fig, axs = plt.subplots(2, 2, figsize=(7.5,7.5))
    fig.suptitle(f"File {file_suffix}",fontsize=20)
    for i,track in enumerate(["cage","dhs","starr","sure"]):
        df_binned=bin_file(df,get_track_num(file_suffix)[track])
        row=int(i/2)
        col=i%2
        sns.histplot(x=f"log10_AF",hue="variant_type", bins=50, data=df_binned,ax=axs[row,col])
        axs[row,col].set_yscale("log")
        # make title smaller
        axs[row,col].set_title(f"Track: {track}",fontsize=10)
        if i in [0,2]:
            axs[row,col].set_ylabel("Count")
        else:
            axs[row,col].set_ylabel("")
        if i in [2,3]:
            axs[row,col].set_xlabel("log10(Allele frequency)")
        else:
            axs[row,col].set_xlabel("")    
    fig.savefig(f"/isdata/alab/people/pcr980/DeepCompare/Effect_size_vs_maf/Plots/effect_size_reverse_distribution_{file_suffix}.pdf")
    plt.close(fig)
    
    


for file_suffix in ["promoters_hepg2","promoters_k562","enhancers_hepg2","enhancers_k562"]:
    plot_file(file_suffix)
    logger.info(f"Plotted {file_suffix}")





# nohup python3 eda_reverse.py > eda_reverse.out &