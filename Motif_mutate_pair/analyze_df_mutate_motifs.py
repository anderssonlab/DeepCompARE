import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# Define a custom color map
colors = ["blue", "lightgrey", "yellow"]  # blue for negative, lightgray for zero, red for positive
n_bins = 100  # Increase this number to make the transition smoother
cmap_name = "custom"
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)



def plot_ism2_vs_ism2_wo_protein1(file_name,file_prefix):
    df=pd.read_csv(file_name)
    df["ism2_wo_protein1"]=df['pred_mut1']-df['pred_mut_both']
    df["ism2_wo_protein1_null"]=df['pred_mut1_null']-df['pred_mut_both_null']
    df["distance"]=df["start_rel2"]-df["start_rel1"]
    df["distance_null"]=df["start_rel2"]-df["start_rel_null"]
    for tf in df["protein2"].unique():
        df_sub=df[df["protein2"]==tf].copy()
        if df_sub.shape[0]>0:
            # TODO: LM
            # plot 2 subplots in two columns
            plt.figure(figsize=(15, 4))
            fig,axs=plt.subplots(1,2)
            # first subplot
            sns.scatterplot(data=df_sub, x='ism2_wo_protein1', y='ism_score_mut2',hue="distance",s=10,palette=custom_cmap,ax=axs[0])
            min_val=min(df_sub['ism2_wo_protein1'].min(),df_sub['ism_score_mut2'].min())
            max_val=max(df_sub['ism2_wo_protein1'].max(),df_sub['ism_score_mut2'].max())
            axs[0].plot([min_val,max_val],[min_val,max_val],color='black',linestyle='--',linewidth=0.5)
            axs[0].set_title(f"{tf},{file_prefix}")
            axs[0].legend(loc='lower right')
            # second subplot: all about null
            sns.scatterplot(data=df_sub, x='ism2_wo_protein1_null', y='ism_score_mut2',hue="distance_null",s=10,palette=custom_cmap,ax=axs[1])
            min_val=min(df_sub['ism2_wo_protein1_null'].min(),df_sub['ism_score_mut2'].min())
            max_val=max(df_sub['ism2_wo_protein1_null'].max(),df_sub['ism_score_mut2'].max())
            axs[1].plot([min_val,max_val],[min_val,max_val],color='black',linestyle='--',linewidth=0.5)
            axs[1].set_title(f"null")
            axs[1].legend(loc='lower right')
            plt.savefig(f"Plots/{tf}_{file_prefix}.pdf",dpi=300)
            plt.close()
            

# TODO: select a distance threshold (or gradient of distance threholds), 
# use wilcoxon signed rank test to prove that ism2_wo_protein1 is the same as ism2 above the distance threshold 


if __name__ == "__main__":
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_enhancers_hepg2_remap_Hep-G2,track4.csv",
                                  "enhancers_hepg2_4")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_enhancers_k562_remap_K-562,track5.csv",
                                  "enhancers_k562_5")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_promoters_hepg2_remap_Hep-G2,track0.csv",
                                  "promoters_hepg2_0")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_promoters_k562_remap_K-562,track1.csv",
                                  "promoters_k562_1")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_promoters_hepg2_remap_Hep-G2,track6.csv",
                                  "promoters_hepg2_6")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_promoters_k562_remap_K-562,track7.csv",
                                  "promoters_k562_7")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_enhancers_hepg2_remap_Hep-G2,track6.csv",
                                  "enhancers_hepg2_6")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_enhancers_k562_remap_K-562,track7.csv",
                                  "enhancers_k562_7")




# nohup python3 analyze_df_mutate_motifs.py &


