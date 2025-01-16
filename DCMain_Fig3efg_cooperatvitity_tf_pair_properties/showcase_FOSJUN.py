import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import change_track_name


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


def read_file(file_name):
    df=pd.read_csv(file_name)
    df["distance"]=np.abs(df["start_rel2"]-df["start_rel1"])
    if "hepg2" in file_name:
        track_list=[0,2,4,6]
    elif "k562" in file_name:
        track_list=[1,3,5,7]
    else:
        raise ValueError("Invalid dataset")
    cols_retain=[f"isa1_track{track_num}" for track_num in track_list]+ \
        [f"isa2_track{track_num}" for track_num in track_list]+ \
        [f"isa_both_track{track_num}" for track_num in track_list]+ \
        ["distance"]
    df=df[["protein1","protein2"]+cols_retain]
    for track_num in track_list:
        df[f"isa2_wo_1_track{track_num}"]=df[f"isa_both_track{track_num}"]-df[f"isa1_track{track_num}"]
    df=change_track_name(df)
    return df



def plot(df, tf, track, title, out_file):
    df_sub = df[df["protein2"] == tf].copy()
    plt.figure(figsize=(2.3, 2))
    # Thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    # Use your previous color scheme
    cmap = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)  # Replace 'cherry' with your specific palette name
    # Create scatter plot with color mapping
    scatter = sns.scatterplot(
        data=df_sub, 
        x=f'isa2_wo_1_{track}', 
        y=f'isa2_{track}', 
        hue="distance", 
        s=2, 
        palette=cmap, 
        legend=False
    )
    # Add a continuous color bar
    norm = plt.Normalize(df_sub["distance"].min(), df_sub["distance"].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.ax.tick_params(labelsize=5)
    cbar.set_label("Distance", fontsize=5)
    # Get min and max values
    min_val = min(df_sub[f'isa2_wo_1_{track}'].min(), df_sub[f'isa2_{track}'].min())
    max_val = max(df_sub[f'isa2_wo_1_{track}'].max(), df_sub[f'isa2_{track}'].max())
    # Identify codependent and redundant pairs
    # df_codependent = df_sub[df_sub[f'isa2_wo_1_{track}'] < df_sub[f'isa2_{track}']]
    # df_redundant = df_sub[df_sub[f'isa2_wo_1_{track}'] > df_sub[f'isa2_{track}']]
    # Annotate number of codependent and redundant pairs
    # plt.text(min_val + 0.1, max_val - 0.1, f"{df_codependent.shape[0]}\n upper-diagonal pairs", fontsize=5)
    # plt.text(min_val + 0.2, min_val + 0.1, f"{df_redundant.shape[0]} \nlower-diagonal pairs", fontsize=5)
    # Add diagonal reference line
    plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=0.1)
    # Labels and formatting
    plt.xlabel("SuRE ISA without partner", fontsize=7)
    plt.ylabel("SuRE ISA with partner", fontsize=7)
    xticks = [0,0.5,1,1.5,2]
    yticks = [0,0.5,1,1.5,2]
    plt.xticks(xticks, fontsize=5)
    plt.yticks(yticks, fontsize=5)
    plt.tight_layout()
    # tight margin
    plt.margins(0.01)
    plt.savefig(f"Plots/{out_file}", dpi=300)
    plt.close()




# find TF very different in promoter v.s. enhancer
df_ci_enhancers=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_enhancer.csv")
df_ci_promoters=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_promoter.csv")
df=pd.merge(df_ci_enhancers,df_ci_promoters,on=["protein2"],suffixes=("_enhancer","_promoter"))
df=df[df["c_sum_enhancer"]>10].reset_index(drop=True)
df=df[df["c_sum_promoter"]>10].reset_index(drop=True)
df["diff"]=df["cooperativity_index_enhancer"]-df["cooperativity_index_promoter"]
# sort by difference
df.sort_values("diff",ascending=False,inplace=True)
df.to_csv("temp.csv",index=False)




# read data
df_enhancers_k562=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_k562.csv")
df_promoters_k562=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_k562.csv")

plot(df_enhancers_k562,"FOS::JUN","sure","Enhancers K562", "FOSJUN_enhancers_k562_sure.pdf")
plot(df_promoters_k562,"FOS::JUN","sure","Promoters K562", "FOSJUN_promoters_k562_sure.pdf")