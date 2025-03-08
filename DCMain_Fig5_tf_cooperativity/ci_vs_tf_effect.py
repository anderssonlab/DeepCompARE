import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics



def violin_ci_vs_tf_effect(effect_col, ylabel, output_file):
    """
    Plots a violin plot for the cooperativity index bins against a given effect column,
    annotates with Mann-Whitney U test p-values, and saves the plot to a file.

    Parameters:
        effect_col (str): Column name for the y-axis continuous variable.
        ylabel (str): Label for the y-axis.
        output_file (str): File path to save the plot.
    """
    plot_violin_with_statistics(
        df=df,
        x_col="tf_type",
        y_col=effect_col,
        x_label="TF type",
        y_label=ylabel,
        title=None,
        output_file=output_file
    )




def scatter_ci_vs_tf_effect(df,modality,cell_line):
    plt.figure(figsize=(2.3,2.3))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(x="cooperativity_index",y=f"avg_isa_{modality}_activity",data=df,hue="tf_type",palette={"intermediate":"gray","codependent":"orangered","redundant":"dodgerblue"},s=5)
    # add pearson correlation and p
    r,p=pearsonr(df["cooperativity_index"],df[f"avg_isa_{modality}_activity"])
    # text at center
    plt.text(0.6,0.6,f"r={r:.2f}\np={p:.2e}",horizontalalignment='center',verticalalignment='center',transform=plt.gca().transAxes,fontsize=5)
    # legend
    plt.legend(title="TF type",fontsize=5, title_fontsize=5)
    # ylim upper: max(df[f"avg_isa_{modality}_activity"])+0.1
    plt.ylim(0,max(df[f"avg_isa_{modality}_activity"])+0.1)
    plt.xticks([0,0.3,0.5,0.7,1],fontsize=5)
    plt.yticks(fontsize=5)
    plt.xlabel("Cooperativity index", fontsize=7)
    plt.ylabel(f"Motif ISA ({modality.upper()})", fontsize=7)
    plt.tight_layout()
    plt.savefig(f"ci_vs_{modality}_scatter_{cell_line}.pdf")
    plt.close()







for cell_line in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_{cell_line}_pe.csv")
    #
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}_pe.txt",header=None).iloc[:,0].tolist()
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}_pe.txt",header=None).iloc[:,0].tolist()
    #
    df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_pe.csv")
    df_coop=df_coop[df_coop["c_sum"]>1].reset_index(drop=True)
    # rename protein2 to protein
    df_coop.rename(columns={"protein2":"protein"},inplace=True)
    # merge with df
    df=df.merge(df_coop[["protein","cooperativity_index"]],on="protein")
    df["tf_type"]="intermediate"
    df.loc[df["protein"].isin(tfs_redundant),"tf_type"]="redundant"
    df.loc[df["protein"].isin(tfs_codependent),"tf_type"]="codependent"
    # convert tf_type to category
    df["tf_type"]=pd.Categorical(df["tf_type"],categories=["redundant","intermediate","codependent"],ordered=True)
    scatter_ci_vs_tf_effect(df,"cage",cell_line)
    scatter_ci_vs_tf_effect(df,"starr",cell_line)
    # violin_ci_vs_tf_effect("avg_isa_starr_activity","Individual TF effect on STARR",f"ci_vs_starr_{cell_line}.pdf")
    # violin_ci_vs_tf_effect("avg_isa_cage_activity","Individual TF effect on CAGE",f"ci_vs_cage_{cell_line}.pdf")

