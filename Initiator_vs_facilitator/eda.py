import pandas as pd
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


def read_file(file_suffix,rare_thresh=1e-3,common_thresh=1e-2):
    df=pd.read_csv(f"maf_with_effect_size_{file_suffix}.csv",header=None,index_col=0)
    df.reset_index(drop=True,inplace=True)
    df.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
    df=bin_and_label(df, "AF", [0,rare_thresh,common_thresh,1])
    df["variant_type"]=["rare" if i==f"0 - {rare_thresh}" else "low" if i==f"{rare_thresh} - {common_thresh}" else "common" for i in df["Bin"]]
    df.drop(columns=["Bin"],inplace=True)
    df=df[df["variant_type"]!="low"].reset_index(drop=True)
    df["variant_type"]=pd.Categorical(df["variant_type"],categories=["rare","common"],ordered=True)
    return df


def plot_file(file_suffix):
    df=read_file(file_suffix)
    fig, axs = plt.subplots(2, 2, figsize=(7.5,7.5))
    fig.suptitle(f"SNV effect size in {file_suffix}",fontsize=20)
    for i,track in enumerate(["cage","dhs","starr","sure"]):
        row=int(i/2)
        col=i%2
        sns.histplot(x=f"track_{get_track_num(file_suffix)[track]}",hue="variant_type", bins=50, data=df,ax=axs[row,col])
        axs[row,col].set_yscale("log")
        # make title smaller
        axs[row,col].set_title(f"Track: {track}",fontsize=10)
        if i in [0,2]:
            axs[row,col].set_ylabel("Count")
        else:
            axs[row,col].set_ylabel("")
        if i in [2,3]:
            axs[row,col].set_xlabel("Predicted effect size")
        else:
            axs[row,col].set_xlabel("")    
        

    fig.savefig(f"/isdata/alab/people/pcr980/DeepCompare/Effect_size_vs_maf/Plots/effect_size_distribution_{file_suffix}.pdf")
    plt.close(fig)
    
    


for file_suffix in ["promoters_hepg2","promoters_k562","enhancers_hepg2","enhancers_k562"]:
    plot_file(file_suffix)
    logger.info(f"Plotted {file_suffix}")



# nohup python3 eda.py > eda.out &