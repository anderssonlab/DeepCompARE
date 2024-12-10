import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer
#-------------------------------
# calculate cell type specificity
#-------------------------------


def read_joint_dispersion():
    df_dispersion_hepg2=pd.read_csv("TFs.dispersionEstimates.hepG2.tab",sep="\t")
    df_dispersion_k562=pd.read_csv("TFs.dispersionEstimates.k562.tab",sep="\t")
    df_dispersion=pd.concat([df_dispersion_hepg2,df_dispersion_k562],axis=0).reset_index(drop=True)
    df_dispersion=df_dispersion.drop_duplicates().reset_index(drop=True)
    return df_dispersion



tfs_redundant = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_merged.txt",header=None)[0].to_list()
tfs_codependent = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_merged.txt",header=None)[0].to_list()
tf_redundant=split_dimer(tfs_redundant)
tfs_codependent=split_dimer(tfs_codependent)
df_dispersion=read_joint_dispersion()
df_dispersion["tf_type"]=np.where(df_dispersion["gene"].isin(tfs_redundant),"redundant",np.where(df_dispersion["gene"].isin(tfs_codependent),"codependent","other"))
ci_codependent = df_dispersion[df_dispersion["tf_type"]=="codependent"]
ci_redundant= df_dispersion[df_dispersion["tf_type"]=="redundant"]
stat, p_value = mannwhitneyu(ci_codependent["gini"],ci_redundant["gini"])
print(f"p-value: {p_value}")

df_dispersion["tf_type"]=pd.Categorical(df_dispersion["tf_type"],categories=["redundant","other","codependent"],ordered=True)




# Assuming df_dispersion is your DataFrame
categories = df_dispersion["tf_type"].unique()
p_values = []

# Calculate pairwise p-values
for cat1, cat2 in combinations(categories, 2):
    group1 = df_dispersion[df_dispersion["tf_type"] == cat1]["gini"]
    group2 = df_dispersion[df_dispersion["tf_type"] == cat2]["gini"]
    stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
    p_values.append({"cat1": cat1, "cat2": cat2, "p_value": p_value})

# Convert to DataFrame for easier handling
p_values_df = pd.DataFrame(p_values)

# Count the number of observations in each category
category_counts = df_dispersion["tf_type"].value_counts()

# Plot boxplot
plt.figure(figsize=(6, 6))
sns.boxplot(data=df_dispersion, x="tf_type", y="gini")
plt.title("Gini coefficient of different TF types")

# Annotate plot with p-values
y_max = df_dispersion["gini"].max() + 0.1  # Adjust for annotation placement
annotation_offset = 0.1  # Increase spacing between annotation lines
line_offset = 0.02       # Adjust vertical height of lines
categories_list = categories.tolist()  # To ensure consistent order

for i, row in p_values_df.iterrows():
    cat1, cat2, p_value = row["cat1"], row["cat2"], row["p_value"]
    x1, x2 = categories_list.index(cat1), categories_list.index(cat2)
    y = y_max + i * annotation_offset  # Increase vertical spacing
    plt.plot([x1, x1, x2, x2], [y, y + line_offset, y + line_offset, y], lw=1.5, color="black")  # Draw annotation lines
    plt.text((x1 + x2) / 2, y + line_offset + 0.01, f"p={p_value:.3e}", ha="center")  # Add p-value

# Add counts to x-axis labels
xticks_labels = [f"{cat}\n(n={category_counts[cat]})" for cat in categories]
plt.xticks(ticks=range(len(categories)), labels=xticks_labels)

# Finalize and save the plot
plt.xlabel("TF Type")
plt.ylabel("Gini Coefficient")
plt.tight_layout()
plt.savefig("cell_type_specificity.pdf")
plt.close()
