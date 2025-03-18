import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


from utils_ppi import read_pooled_found_tf_pairs

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42




suffix="pe"

df_bait_pooled=read_pooled_found_tf_pairs()


for bait in ["BACH1", "RFX5", "IKZF1", "MAFG", "RREB1"]:
    for mode in ["linearity_index", "cooperativity_index"]:
        # read cooperativity index
        df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_{suffix}.csv")
        df_coop = assign_cooperativity(df_coop,1,0.9,0.3,0.7)
        
        # read published PPI
        df_ppi = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/DCMain_Fig4_tf_pair_cooperativity/Pd1_Petra_data/2025-03-07_s10_PublishedPPIandProtComplexes_k562_pe.txt", sep='\t')
        df_ppi.rename(columns={'SWI/SNF': 'SWI_SNF'}, inplace=True)
        df_coop = pd.merge(df_coop, df_ppi, on=["protein1", "protein2"], how="left")
        # change Reported_PPI to 0,1: 0 for no, 1 for es
        df_coop["Reported_PPI"] = df_coop["Reported_PPI"].fillna(0)
        df_coop["Reported_PPI"] = df_coop["Reported_PPI"].map({"No": 0, "Yes": 1})
        df_coop = df_coop.dropna(subset=["Reported_PPI"]).reset_index(drop=True)
        # subset to the bait
        df_coop = df_coop[df_coop["protein2"] == bait].reset_index(drop=True)
        
        # get the found TFs in df_bait
        df_bait = df_bait_pooled[df_bait_pooled["protein2"] == bait].reset_index(drop=True)
        df_coop["found"]=df_coop["protein1"].isin(df_bait["protein1"])
        df_sub=df_coop[df_coop["found"]]
        df_sub.to_csv("eda.csv", index=False)
        plt.figure(figsize=(2.5, 2.5)) 
        # thin frame
        plt.gca().spines['top'].set_linewidth(0.5)
        plt.gca().spines['right'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)
        # violin, no color
        ax = sns.violinplot(x="Reported_PPI", y=mode, data=df_coop, inner="quartile", linewidth=0.5,cut=0)

        # Remove fill color from each violin
        for violin in ax.collections:
            violin.set_facecolor('none')  # Remove fill
            violin.set_edgecolor('black')     
            
        # Annotate with protein names for the found candidates
        texts = []
        for _, row in df_coop[df_coop["found"]].iterrows():
            texts.append(plt.text(row["Reported_PPI"], row[mode], row["protein1"], fontsize=5,color='black'))
            
        plt.title(bait, fontsize=7)
        plt.xlabel("Reported PPI", fontsize=7)
        plt.ylabel(mode, fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"eda_{mode}_{bait}_{suffix}.pdf")
        plt.close()
    
        
# coop index v.s. log fold change