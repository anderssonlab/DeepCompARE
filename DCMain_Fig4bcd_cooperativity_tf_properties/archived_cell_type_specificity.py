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



cell_line="k562"

df_dispersion=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/gtex.dispersionEstimates.tab",sep="\t")

df_tf=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}.csv")
# merge with df_dispersion
df_tf=df_tf.merge(df_dispersion,left_on="protein2",right_on="symbol",how="inner")
df_tf.drop_duplicates(subset="protein2",inplace=True)


tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt",header=None).iloc[:,0].tolist()
tfs_codependent=split_dimer(tfs_codependent)
tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt",header=None).iloc[:,0].tolist()
tfs_redundant=split_dimer(tfs_redundant)
df_tf["tf_type"]="other"
df_tf.loc[df_tf["protein2"].isin(tfs_codependent),"tf_type"]="codependent"
df_tf.loc[df_tf["protein2"].isin(tfs_redundant),"tf_type"]="redundant"



# pearson correlation
from scipy.stats import pearsonr
pearsonr(df_tf["adjusted_dispersion"],df_tf["cooperativity_index"])

ci_codependent = df_tf[df_tf["tf_type"]=="codependent"]
ci_redundant= df_tf[df_tf["tf_type"]=="redundant"]
stat, p_value = mannwhitneyu(ci_codependent["adjusted_dispersion"],ci_redundant["adjusted_dispersion"])
print(f"p-value: {p_value}")

df_tf["tf_type"]=pd.Categorical(df_tf["tf_type"],categories=["redundant","other","codependent"],ordered=True)




categories = df_tf["tf_type"].cat.categories
p_values = []

# Calculate pairwise p-values
for cat1, cat2 in combinations(categories, 2):
    group1 = df_tf[df_tf["tf_type"] == cat1]["gini"]
    group2 = df_tf[df_tf["tf_type"] == cat2]["gini"]
    stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
    p_values.append({"cat1": cat1, "cat2": cat2, "p_value": p_value})

# Convert to DataFrame for easier handling
p_values_df = pd.DataFrame(p_values)

# Count the number of observations in each category
category_counts = df_tf["tf_type"].value_counts()

# Plot boxplot
plt.figure(figsize=(6, 6))
sns.boxplot(data=df_tf, x="tf_type", y="gini")
plt.title("Gini coefficient of different TF types")

# Annotate plot with p-values
y_max = df_tf["gini"].max() + 0.1  # Adjust for annotation placement
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
