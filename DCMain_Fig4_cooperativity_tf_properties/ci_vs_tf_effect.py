import pandas as pd



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label
from plotting import plot_violin_with_statistics

cell_line="k562"

df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}.csv")
df_coop=df_coop[df_coop["c_sum"]>1].reset_index(drop=True)
df_effect=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_{cell_line}.csv")


df=pd.merge(df_coop,df_effect,left_on="protein2",right_on="protein",how="inner")
df=bin_and_label(df,"cooperativity_index",[0,0.3,0.7,1])


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
        x_col="cooperativity_index_bin",
        y_col=effect_col,
        x_label="Cooperativity index bin",
        y_label=ylabel,
        title="Cooperativity Index vs. TF Effect",
        output_file=output_file
    )



# plot_ci_vs_tf_effect("avg_isa_starr_activity","Individual TF effect on STARR",f"ci_vs_starr_{cell_line}.pdf")
# plot_ci_vs_tf_effect("avg_isa_cage_activity","Individual TF effect on CAGE",f"ci_vs_cage_{cell_line}.pdf")

plot_ci_vs_tf_effect("avg_isa_dhs_activity","Individual TF effect on DHS",f"ci_vs_dhs_{cell_line}.pdf")
plot_ci_vs_tf_effect("avg_isa_sure_activity","Individual TF effect on sure",f"ci_vs_sure_{cell_line}.pdf")



