import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# adjust text


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity



for bait in ["BACH1", "RFX5", "IKZF1", "MAFG", "RREB1"]:
    for mode in ["linearity_index", "cooperativity_index"]:
        # preprocess
        df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_pe.csv")
        df_coop = assign_cooperativity(df_coop,0.3,0.7)

        df_ppi = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/DCMain_Fig4_tf_pair_cooperativity/Pd1_Petra_data/2025-03-07_s10_PublishedPPIandProtComplexes_k562_pe.txt", sep='\t')
        df_ppi.rename(columns={'SWI/SNF': 'SWI_SNF'}, inplace=True)
        df_coop = pd.merge(df_coop, df_ppi, on=["protein1", "protein2"], how="left")
        
        # change Reported_PPI to 0,1: 0 for no, 1 for yes
        df_coop["Reported_PPI"] = df_coop["Reported_PPI"].fillna(0)
        df_coop["Reported_PPI"] = df_coop["Reported_PPI"].map({"No": 0, "Yes": 1})
        df_coop = df_coop.dropna(subset=["Reported_PPI"]).reset_index(drop=True)
        
        df_coop_bait = df_coop[df_coop["protein2"] == bait].reset_index(drop=True)
        df_exp = pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
        
        # Which genes in df_exp are in df_coop_bait?
        candidates = df_coop_bait["protein1"].tolist()
        # Split dimers in candidates
        candidates = [x.split("::") for x in candidates]
        # Flatten the list
        candidates = [item for sublist in candidates for item in sublist]
        tfs_found = set(df_exp["gene"].values).intersection(candidates)
        df_coop_bait_found = df_coop_bait[df_coop_bait["protein1"].isin(tfs_found)]
        
        plt.figure(figsize=(2.5, 2.5)) 
        # thin frame
        plt.gca().spines['top'].set_linewidth(0.5)
        plt.gca().spines['right'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)
        # violin
        sns.violinplot(x="Reported_PPI", y=mode, data=df_coop_bait, inner="quartile",linewidth=0.5)
    
        # Annotate with protein names for the found candidates
        texts = []
        for _, row in df_coop_bait_found.iterrows():
            texts.append(plt.text(row["Reported_PPI"], row[mode], row["protein1"], fontsize=5,color='black'))
            
        plt.title(bait, fontsize=7)
        plt.xlabel("Reported PPI", fontsize=7)
        plt.ylabel(mode, fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"eda_{mode}_{bait}.pdf")
        plt.close()
    
        
# coop index v.s. log fold change