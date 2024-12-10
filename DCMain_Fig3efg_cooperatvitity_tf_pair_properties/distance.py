import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import pearsonr,mannwhitneyu
from adjustText import adjust_text

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label

cell_line="merged"

df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}.csv")
df=df[df["c_sum"]>1].reset_index(drop=True)
pearsonr(df["cooperativity_index"],df["distance"])[0]

df=bin_and_label(df,"cooperativity_index",[0,0.3,0.7,1])


#-------------------------
# Analysis1: plot distance distribution of cooperativity index bin
#-------------------------
# Perform pairwise Mann-Whitney U tests
# bins are the catepories of cooperativity index
bins = df["cooperativity_index_bin"].cat.categories.astype(str)
bin_counts = df["cooperativity_index_bin"].value_counts()

# Perform pairwise Mann-Whitney U tests
p_values = []
for bin1, bin2 in combinations(bins, 2):
    group1 = df[df["cooperativity_index_bin"] == bin1]["distance"]
    group2 = df[df["cooperativity_index_bin"] == bin2]["distance"]
    stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
    p_values.append({"bin1": bin1, "bin2": bin2, "p_value": p_value})

p_values_df = pd.DataFrame(p_values)
significant_results = p_values_df[p_values_df["p_value"] < 0.05]

# Plot violin plot
plt.figure(figsize=(6, 6))
sns.violinplot(data=df, x="cooperativity_index_bin", y="distance", cut=0)

# Annotate plot with p-values
for i, row in significant_results.iterrows():
    bin1, bin2 = row["bin1"], row["bin2"]
    p_value = row["p_value"]
    y_max = max(df["distance"]) + 10  # Adjust for annotation placement
    x1, x2 = bins.tolist().index(bin1), bins.tolist().index(bin2)
    y = y_max + i * 20  # Adjust for stacking annotations
    plt.plot([x1, x1, x2, x2], [y, y + 1, y + 1, y], lw=1.5, color="black")
    plt.text((x1 + x2) / 2, y + 4, f"p={p_value:.3e}", ha="center")

# Add counts to x-axis labels
bin_labels = [f"{bin_}\n(n={bin_counts[bin_]})" for bin_ in bins]
plt.xticks(ticks=range(len(bins)), labels=bin_labels)

# Finalize the plot
plt.xlabel("Cooperativity index bin")
plt.ylabel("Median distance between TF pair (bp)")
plt.tight_layout()
plt.savefig(f"distance_vs_cooperativity_index_bin_{cell_line}.pdf")
plt.close()

#-------------------------
# Analysis 2: compare distance between codependent and redundant pairs
#-------------------------
df["cooperativity"]=np.where(df["cooperativity_index"]>0.7,"codependent",np.where(df["cooperativity_index"]<0.3,"redundant",""))


test_results = []
for protein in df["protein2"].unique():
    # Subset data for the current protein
    df_subset = df[df["protein2"] == protein].reset_index(drop=True)
    redundant_distances = df_subset.loc[df_subset["cooperativity"] == "redundant", "distance"]
    codependent_distances = df_subset.loc[df_subset["cooperativity"] == "codependent", "distance"]
    if len(redundant_distances) >= 2 and len(codependent_distances) >= 2:
        stat, p_value = mannwhitneyu(redundant_distances, codependent_distances, alternative="two-sided")
        test_results.append({"protein": protein, 
                             "statistic": stat, 
                             "p_value": p_value,
                             "redundant_median": np.median(redundant_distances),
                             "codependent_median": np.median(codependent_distances)})

# Convert test results into a DataFrame
df_res = pd.DataFrame(test_results)
df_res["significant"] = df_res["p_value"] < 0.05
df_res["significant"].sum()
# scatter plot
plt.figure(figsize=(6,6))
sns.scatterplot(data=df_res,x="redundant_median",y="codependent_median",hue="significant")
plt.xlabel("Median distance of redundant pairs")
plt.ylabel("Median distance of codependent pairs")
# add diagonal line
min_val=min(df_res[["codependent_median","redundant_median"]].min())
max_val=max(df_res[["codependent_median","redundant_median"]].max())
plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
# annotate each tf
texts = []
for i in range(df_res.shape[0]):
    if np.abs([df_res["redundant_median"][i]-df_res["codependent_median"][i]])>60:
        texts.append(plt.text(df_res["redundant_median"][i],df_res["codependent_median"][i],df_res["protein"][i]))

adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))

plt.savefig(f"distance_codependent_vs_redundant_{cell_line}.pdf")
plt.close()