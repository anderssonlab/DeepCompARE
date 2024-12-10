import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label



# TODO: use "_merged" data, use motif_info_400
df_coop=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_k562.csv")
df_coop=df_coop[df_coop["c_sum"]>1].reset_index(drop=True)
df_effect=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd9_TF_effect_and_constraint/tf_effect_and_constraints_k562.csv")


df=pd.merge(df_coop,df_effect,left_on="protein2",right_on="protein",how="inner")
df=bin_and_label(df,"cooperativity_index",[0,0.3,0.7,1])


#-------------------------
# Analysis1: plot distance distribution of cooperativity index bin
#-------------------------
# Perform pairwise Mann-Whitney U tests
# bins are the catepories of cooperativity index

def plot_ci_vs_tf_effect(effect_col,ylabel,output_file):
    bins = df["cooperativity_index_bin"].cat.categories.astype(str)
    bin_counts = df["cooperativity_index_bin"].value_counts()
    # Perform pairwise Mann-Whitney U tests
    p_values = []
    for bin1, bin2 in combinations(bins, 2):
        group1 = df[df["cooperativity_index_bin"] == bin1][effect_col]
        group2 = df[df["cooperativity_index_bin"] == bin2][effect_col]
        stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
        p_values.append({"bin1": bin1, "bin2": bin2, "p_value": p_value})
    p_values_df = pd.DataFrame(p_values)
    # Plot violin plot
    plt.figure(figsize=(6, 6))
    sns.violinplot(data=df, x="cooperativity_index_bin", y=effect_col, cut=0)
    # Annotate plot with p-values
    for i, row in p_values_df.iterrows():
        bin1, bin2 = row["bin1"], row["bin2"]
        p_value = row["p_value"]
        y_max = max(df[effect_col]) + 0.1  # Adjust for annotation placement
        x1, x2 = bins.tolist().index(bin1), bins.tolist().index(bin2)
        y = y_max + i * 0.1  # Adjust for stacking annotations
        plt.plot([x1, x1, x2, x2], [y, y + 0.05, y + 0.05, y], lw=1.5, color="black")
        plt.text((x1 + x2) / 2, y + 0.06, f"p={p_value:.3e}", ha="center")
    # Add counts to x-axis labels
    bin_labels = [f"{bin_}\n(n={bin_counts[bin_]})" for bin_ in bins]
    plt.xticks(ticks=range(len(bins)), labels=bin_labels)
    # Finalize the plot
    plt.xlabel("Cooperativity index bin")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()





plot_ci_vs_tf_effect("avg_ism_starr_activity","Individual TF effect on STARR","ci_vs_starr.png")
plot_ci_vs_tf_effect("avg_ism_cage_activity","Individual TF effect on CAGE","ci_vs_cage.png")




