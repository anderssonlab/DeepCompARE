
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label


cell_line="k562"

df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}.csv")
df=df[df['c_sum']>1].reset_index(drop=True)

df_chip_coloc=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_{cell_line}_summary.csv")
df_chip_coloc.rename(columns={'tf1':'protein1','tf2':'protein2'},inplace=True)

df=pd.merge(df,df_chip_coloc,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')
df.rename(columns={'mean_abs':'mean_chip_peak_distance',
                   'tf_pair_count':'chip_overlap_count',
                   },inplace=True)

df=bin_and_label(df,"cooperativity_index",[0,0.3,0.7,1])


# ----------------- Plotting -----------------
# Plot chip peak distance vs cooperativity index


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
plt.ylabel("Mean distance between ChIP peaks (bp)")
plt.tight_layout()
plt.savefig(f"chip_distance_vs_ci_{cell_line}.pdf")
plt.close()
