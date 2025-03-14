import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from adjustText import adjust_text


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer
from tf_cooperativity import assign_cooperativity


# Use a fixed colormap
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cmap = mpl.cm.coolwarm





bait_list = ["BACH1", "IKZF1", "MAFG", "RFX5", "RREB1"]

for bait in bait_list:
    print(bait)
    df_bait = pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
    df_htfs = pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv", index_col=0)
    df_htfs=df_htfs[df_htfs["Is TF?"]=="Yes"].reset_index(drop=True)
    # select K562 expressed TFs
    proteins_expressed=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv",sep='\t',header=None).iloc[:,0].values
    df_htfs=df_htfs[df_htfs["HGNC symbol"].isin(proteins_expressed)].reset_index(drop=True)
    # read cooperativity
    df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_pe.csv")
    df_coop = assign_cooperativity(df_coop, 0.3, 0.7)
    # select baits
    df_coop = df_coop[df_coop["protein2"] == bait].reset_index(drop=True)
    # determine if TF is linear or nonlinear
    df_nonlinear = df_coop[df_coop["cooperativity"]!="Linear"].reset_index(drop=True)
    df_linear = df_coop[df_coop["cooperativity"]=="Linear"].reset_index(drop=True)
    df_bait["is_tf"] = df_bait["gene"].isin(df_htfs["HGNC symbol"])
    df_bait["is_nonlinear"] = df_bait["gene"].isin(df_nonlinear["protein1"])
    df_bait["is_linear"] = df_bait["gene"].isin(df_linear["protein1"])
    df_bait["is_not_investigated"] = ~df_bait["gene"].isin(df_coop["protein1"])
    df_bait["-log10_pvalue"] = df_bait["pvalue"].apply(lambda x: -1 if x == 0 else -1 * np.log10(x))

    
    plt.figure(figsize=(2.5,2.5))
    # Adjust thin frame
    for spine in plt.gca().spines.values():
        spine.set_linewidth(0.5)
    
    # Plot everything as background
    sns.scatterplot(
        data=df_bait,
        x="logFC",
        y="-log10_pvalue",
        s=5,
        color="LightGray",
        legend=False,
        alpha=0.3
    )
    
    # Plot uninvestigated TFs
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
    
    # plot lineari TFs 
    df_linear=pd.merge(df_linear,df_bait,how="inner",left_on="protein1",right_on="gene")
    sns.scatterplot(data=df_linear,x="logFC",y="-log10_pvalue",s=5,color="black",label="Linear TFs",marker="x")
    
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

        
    # Add gene names for investigated TFs
    texts = []
    for i in range(len(df_nonlinear)):
        texts.append(plt.text(df_nonlinear.loc[i, "logFC"],
                 df_nonlinear.loc[i, "-log10_pvalue"],
                 df_nonlinear.loc[i, "gene"],
                 fontsize=5))
    adjust_text(texts)
    
    # Add horizontal significance line
    plt.axhline(-1 * np.log10(0.05), color="black", linestyle="--", linewidth=0.5)
    plt.xlabel("log fold change", fontsize=7)
    plt.ylabel("-log10 pvalue", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    # legend font size: 5
    plt.legend(fontsize=5)
    plt.title(f"{bait}", fontsize=7)
    # Add the same coolwarm color bar to each plot
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm)
    cbar.set_label("Cooperativity Index", fontsize=5)  # Adjust colorbar title size
    cbar.ax.tick_params(labelsize=5)  # Adjust colorbar tick size 
    plt.tight_layout()
    plt.savefig(f"volcano_{bait}.pdf")
    plt.close()


