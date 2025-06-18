import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from adjustText import adjust_text
import matplotlib.lines as mlines



from utils_ppi import get_htfs_list


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity





mpl.rcParams['pdf.fonttype']=42



# Use a fixed colormap
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cmap = mpl.cm.coolwarm





bait_list = ["BACH1", "IKZF1", "MAFG", "RFX5", "RREB1"]


for bait in bait_list:
    df_bait = pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
    htfs=get_htfs_list()
    # select K562 expressed TFs
    proteins_expressed=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv",sep='\t',header=None).iloc[:,0].values
    htfs=list(set(htfs).intersection(proteins_expressed))
    # read cooperativity
    df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_pe.csv")
    df_coop = assign_cooperativity(df_coop, 1, 0.9, 0.3, 0.7)
    # select baits
    df_coop = df_coop[df_coop["protein2"] == bait].reset_index(drop=True)
    # determine if TF is linear or nonlinear
    df_nonlinear = df_coop[df_coop["cooperativity"]!="Independent"].reset_index(drop=True)
    df_linear = df_coop[df_coop["cooperativity"]=="Independent"].reset_index(drop=True)
    df_bait["is_tf"] = df_bait["gene"].isin(htfs)
    df_bait["is_nonlinear"] = df_bait["gene"].isin(df_nonlinear["protein1"])
    df_bait["is_linear"] = df_bait["gene"].isin(df_linear["protein1"])
    df_bait["is_not_investigated"] = ~df_bait["gene"].isin(df_coop["protein1"])
    df_bait["-log10_pvalue"] = df_bait["pvalue"].apply(lambda x: -1 if x == 0 else -1 * np.log10(x))

    
    plt.figure(figsize=(1.9,1.9)) # (2.8,2.6) for main. (1.9,1.9) for sup
    # Adjust thin frame
    for spine in plt.gca().spines.values():
        spine.set_linewidth(0.5)
    
    # Plot all tfs as light gray dots
    sns.scatterplot(
        data=df_bait,
        x="logFC",
        y="-log10_pvalue",
        s=5,
        color="LightGray",
        legend=False,
        alpha=0.3
    )
    
    # Plot uninvestigated TFs as black dots
    df_bait_tfs = df_bait[(df_bait["is_tf"]) & (df_bait["is_not_investigated"])].reset_index(drop=True)
    sns.scatterplot(
        data=df_bait_tfs,
        x="logFC",
        y="-log10_pvalue",
        s=5,
        color="black",
        legend=False,
        alpha=0.3
    )
    
    # plot linear TFs as black x
    df_linear=pd.merge(df_linear,df_bait,how="inner",left_on="protein1",right_on="gene")
    sns.scatterplot(data=df_linear,x="logFC",y="-log10_pvalue",s=5,color="black",label="Independent TF partner",marker="x")
    
    # Plot nonlinear TFs with hue based on cooperativity_index
    df_nonlinear=pd.merge(df_nonlinear,df_bait,how="inner",left_on="protein1",right_on="gene")
    # Create scatter plot
    sns.scatterplot(
        data=df_nonlinear,
        x="logFC",
        y="-log10_pvalue",
        hue="cooperativity_index",
        palette=cmap,  # Use the fixed colormap
        hue_norm=norm,  # Apply fixed normalization
        s=5,
        edgecolor=None,
        legend=False
    )

        
    # Add gene names for investigated TFs (only for main)
    # texts = []
    # for i in range(len(df_nonlinear)):
    #     texts.append(plt.text(df_nonlinear.loc[i, "logFC"],
    #              df_nonlinear.loc[i, "-log10_pvalue"],
    #              df_nonlinear.loc[i, "gene"],
    #              fontsize=5))
    # adjust_text(texts)
    
    # Add horizontal significance line
    plt.axhline(-1 * np.log10(0.05), color="black", linestyle="--", linewidth=0.5)
    plt.xlabel("Log fold change", fontsize=7)
    plt.ylabel("-Log10 pvalue", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    # legend font size: 5
    # Custom legend
    legend_elements = [
        # Non-TFs: light gray solid dot
        mlines.Line2D([], [], color="LightGray", marker='o', linestyle='None', markersize=5, label='Non-TF interactors', alpha=0.3),
        # Uninvestigated TFs: black solid dot
        mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=5, label='Uninvestigated TF', alpha=0.3),
        # Independent TF partners: black x
        mlines.Line2D([], [], color='black', marker='x', linestyle='None', markersize=5, label='Independent TF', alpha=1),
        # Cooperative TF partner: hollow circle (uses colormap in plot)
        mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=5, alpha=1, markerfacecolor='none', label='Cooperative TF'),
    ]


    plt.legend(handles=legend_elements, fontsize=5, loc='upper left')

    plt.title(f"{bait} interactors", fontsize=7)
    # # Add the same coolwarm color bar to each plot (Only in main)
    # sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    # cbar = plt.colorbar(sm)
    # cbar.set_label("Synergy score", fontsize=5)  # Adjust colorbar title size
    # cbar.ax.tick_params(labelsize=5)  # Adjust colorbar tick size 
    plt.tight_layout()
    plt.savefig(f"volcano_{bait}.pdf")
    plt.close()


