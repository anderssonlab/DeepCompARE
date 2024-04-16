import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



def plot_ism2_vs_ism2_wo_protein1(file_name,file_prefix):
    df=pd.read_csv(file_name)
    df["ism2_wo_protein1"]=df['pred_mut1']-df['pred_mut_both']
    df["distance"]=df["start_rel2"]-df["start_rel1"]
    for tf in df["protein2"].unique():
        df_sub=df[df["protein2"]==tf].copy()
        if df_sub.shape[0]>0:
            # TODO: GLM
            plt.figure(figsize=(5, 4))
            sns.scatterplot(data=df_sub, x='ism2_wo_protein1', y='ism_score_mut2',hue="distance",s=10)
            min_val=min(df_sub['ism2_wo_protein1'].min(),df_sub['ism_score_mut2'].min())
            max_val=max(df_sub['ism2_wo_protein1'].max(),df_sub['ism_score_mut2'].max())
            plt.plot([min_val,max_val],[min_val,max_val],color='black',linestyle='--',linewidth=0.5)
            plt.title(f"{tf},{file_prefix}")
            plt.legend(loc='lower right')
            plt.savefig(f"Plots/{tf}_{file_prefix}.pdf",dpi=300)
            plt.close()


if __name__ == "__main__":
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_enhancers_hepg2_remap_Hep-G2,track4.csv",
                                  "enhancers_hepg2")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_enhancers_k562_remap_K-562,track5.csv",
                                  "enhancers_k562")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_promoters_hepg2_remap_Hep-G2,track0.csv",
                                  "promoters_hepg2")
    plot_ism2_vs_ism2_wo_protein1("df_mutate_pair_promoters_k562_remap_K-562,track1.csv",
                                  "promoters_k562")




# nohup python3 analyze_df_mutate_motifs.py &


