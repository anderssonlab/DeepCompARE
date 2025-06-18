import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity







import matplotlib
matplotlib.rcParams['pdf.fonttype']=42





def scatter_ci_vs_tf_effect(df,modality,cell_line,suffix):
    df_nonlinear=df[df["cooperativity"]!="Independent"].reset_index(drop=True)
    df_nonlinear["cooperativity"]=pd.Categorical(df_nonlinear["cooperativity"],categories=["Redundant","Intermediate","Synergistic"],ordered=True)

    fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [5, 1]},sharey=True,figsize=(2.3,1.6),constrained_layout=True)
    for ax in axes:  # Loop through both subplots
        ax.spines['top'].set_linewidth(0.5)
        ax.spines['right'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)

    sns.scatterplot(x="cooperativity_index",
                    y=f"avg_isa_{modality}_activity",
                    data=df_nonlinear,
                    hue="cooperativity",
                    ax=axes[0],
                    palette={"Intermediate":"gray","Synergistic":"#d62728","Redundant":"#1f77b4"},s=5)
    # add pearson correlation and p
    r,p=pearsonr(df_nonlinear["cooperativity_index"],df_nonlinear[f"avg_isa_{modality}_activity"])
    # text at center
    axes[0].text(0.5,0.4,f"Pearson R={r:.2f}\nP={p:.2e}",horizontalalignment='center',verticalalignment='center',transform=axes[0].transAxes,fontsize=5)
    
    # legend
    handles, labels = axes[0].get_legend_handles_labels()
    # Create a solid black dot (marker) for Independent TFs
    independent_dot = Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=5, label='Independent')
    # Add it to the legend
    axes[0].legend(handles + [independent_dot], labels + ['Independent'],
               title="TF type", fontsize=5, title_fontsize=5)
    
    axes[0].set_xticks([0, threshold_dict[cell_line][suffix][0], threshold_dict[cell_line][suffix][1], 1])
    axes[0].set_xticklabels([0, threshold_dict[cell_line][suffix][0], threshold_dict[cell_line][suffix][1], 1], fontsize=5)
    axes[0].set_yticks(np.arange(0, 1.2, 0.2))
    axes[0].set_yticklabels([f"{x:.1f}" for x in np.arange(0, 1.2, 0.2)], fontsize=5)
    axes[0].set_xlabel("Synergy score", fontsize=7)
    axes[0].set_ylabel(f"Motif ISA ({modality.upper()})", fontsize=7)
    # strip plot for linear TFs
    sns.stripplot(
        x="cooperativity",
        y=f"avg_isa_{modality}_activity",
        data=df[df["cooperativity"] == "Independent"].reset_index(drop=True),
        ax=axes[1],
        color="black",
        size=2)

    axes[1].set_xticks("")  # Remove x-ticks
    axes[1].set_xlabel("")  # Remove x-label
    ymax = df[f"avg_isa_{modality}_activity"].max() * 1.3
    ymin = df[f"avg_isa_{modality}_activity"].min() * 1.2
    axes[0].set_ylim(ymin, ymax)
    axes[1].set_ylim(ymin, ymax)

    plt.savefig(f"ci_vs_{modality}_scatter_{cell_line}_{suffix}.pdf",dpi=300)
    plt.close()















threshold_dict={"hepg2":{"pe":[0.3,0.7],"dhs":[0.48,0.78]},
                "k562":{"pe":[0.3,0.7],"dhs":[0.44,0.81]}}


for cell_line in ["hepg2","k562"]:
    for suffix in ["pe"]:
        df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_{cell_line}_{suffix}.csv")
        df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_{suffix}.csv")
        df_coop=assign_cooperativity(df_coop,5,0.95,threshold_dict[cell_line][suffix][0],threshold_dict[cell_line][suffix][1])
        # rename protein2 to protein
        df_coop.rename(columns={"protein2":"protein"},inplace=True)
        # merge with df
        df=df.merge(df_coop,on="protein")
        scatter_ci_vs_tf_effect(df,"cage",cell_line,suffix)
        scatter_ci_vs_tf_effect(df,"starr",cell_line,suffix)
