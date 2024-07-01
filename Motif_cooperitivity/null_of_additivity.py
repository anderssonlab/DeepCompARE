""" Understand how TF cooperitivity die down as distance increases"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_null(track_name,hepg2_track,k562_track):
    
    df_enhancers_hepg2=pd.read_csv(f"df_mutate_pair_enhancers_hepg2_remap_Hep-G2,track{hepg2_track}.csv")
    df_enhancers_k562=pd.read_csv(f"df_mutate_pair_enhancers_k562_remap_K-562,track{k562_track}.csv")
    df_promoters_hepg2=pd.read_csv(f"df_mutate_pair_promoters_hepg2_remap_Hep-G2,track{hepg2_track}.csv")
    df_promoters_k562=pd.read_csv(f"df_mutate_pair_promoters_k562_remap_K-562,track{k562_track}.csv")

    df_enhancers_hepg2["data"]="enhancers_hepg2"
    df_enhancers_k562["data"]="enhancers_k562"
    df_promoters_hepg2["data"]="promoters_hepg2"
    df_promoters_k562["data"]="promoters_k562"

    df=pd.concat([df_enhancers_hepg2,
                df_enhancers_k562,
                df_promoters_hepg2,
                df_promoters_k562])

    df["ism2_wo_protein1"]=df['pred_mut1']-df['pred_mut_both']
    df["diff"]=np.abs(df["ism_score_mut2"]-df["ism2_wo_protein1"])
    df["ism2_wo_protein1_null"]=df['pred_mut1_null']-df['pred_mut_both_null']
    df["diff_null"]=np.abs(df["ism_score_mut2"]-df["ism2_wo_protein1_null"])
    df["distance"]=np.abs(df["start_rel2"]-df["start_rel1"])
    df["distance_null"]=np.abs(df["start_rel2"]-df["start_rel_null"])

    
    # group df by distance and data, calculate mean of diff and diff_null
    
    df_alt=df.groupby(["distance","data"]).agg({"diff":"mean"}).reset_index()
    df_null=df.groupby(["distance_null","data"]).agg({"diff_null":"mean"}).reset_index()
    
    
    # plot line plot using seaborn for distance vs diff
    sns.lineplot(x="distance", y="diff", data=df_alt, hue="data")
    sns.lineplot(x="distance_null", y="diff_null", data=df_null, hue="data", alpha=0.3,legend=False)
    plt.title(track_name)
    plt.savefig(f"distance_vs_diff_{track_name}.png")
    plt.close()
            
plot_null("cage_track",0,1)
plot_null("dhs_track",2,3)
plot_null("starr_track",4,5)
plot_null("sure_track",6,7)

# nohup python3 null_hypothesis.py > null_hypothesis.out &

