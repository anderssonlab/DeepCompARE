import pandas as pd



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics

cell_line="k562"

df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_{cell_line}.csv")

tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt",header=None).iloc[:,0].tolist()
tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt",header=None).iloc[:,0].tolist()

df["tf_type"]="other"
df.loc[df["protein"].isin(tfs_redundant),"tf_type"]="redundant"
df.loc[df["protein"].isin(tfs_codependent),"tf_type"]="codependent"


# convert tf_type to category
df["tf_type"]=pd.Categorical(df["tf_type"],categories=["redundant","other","codependent"],ordered=True)
#-------------------------
# Analysis1: plot distance distribution of cooperativity index bin
#-------------------------

def plot_ci_vs_tf_effect(effect_col, ylabel, output_file):
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



plot_ci_vs_tf_effect("avg_isa_starr_activity","Individual TF effect on STARR",f"ci_vs_starr_{cell_line}.pdf")
plot_ci_vs_tf_effect("avg_isa_cage_activity","Individual TF effect on CAGE",f"ci_vs_cage_{cell_line}.pdf")

