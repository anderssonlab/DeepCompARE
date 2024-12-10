import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations



#---------------------------------
# Read and merge data
#---------------------------------

df_coop=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_merged.csv")
df_coop=df_coop[df_coop["c_sum"]>1].reset_index(drop=True)

df_ppi=pd.read_csv("/isdata/alab/people/pcr980/Resource/2024-10-24_s1_PublishedPPIannotation_withProtComplexes.txt",sep='\t')
df_ppi.rename(columns={'SWI/SNF_lowerTh':'SWI_SNF_lowerTh',
                       'SWI/SNF':'SWI_SNF',
                       },inplace=True)
# remove columns protPair','coopIndex_avg', 'sum_ci_avg'
df_ppi.drop(columns=['protPair','coopIndex_avg', 'sum_ci_avg'],inplace=True)
df_ppi2=df_ppi.copy()
df_ppi2.rename(columns={'protein1':'protein2','protein2':'protein1'},inplace=True)
df_ppi=pd.concat([df_ppi,df_ppi2],axis=0).reset_index(drop=True)

df=pd.merge(df_coop,df_ppi,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')



#---------------------------------
# 1. Plot distribution of CI split by Reported_PPI
#---------------------------------

# Perform Mann-Whitney U test
yes_values = df[df['Reported_PPI'] == 'Yes']['cooperativity_index']
no_values = df[df['Reported_PPI'] == 'No']['cooperativity_index']
stat, p_value = mannwhitneyu(yes_values, no_values, alternative='two-sided')
counts = df['Reported_PPI'].value_counts()
# plot
plt.figure(figsize=(4, 4))
sns.violinplot(data=df, x='Reported_PPI', y='cooperativity_index', cut=0)
# Add the p-value annotation
y_max = df['cooperativity_index'].max() + 0.2  # Adjust for annotation placement
x1, x2 = 0, 1  # Positions of the two categories on the x-axis
plt.plot([x1, x1, x2, x2], [y_max, y_max + 0.02, y_max + 0.02, y_max], lw=1.5, color="black")  # Line
plt.text((x1 + x2) / 2, y_max + 0.03, f"p={p_value:.3e}", ha='center')
# Update x-axis tick labels to include counts
plt.xticks(ticks=[0, 1], labels=[f"No\n(n={counts['No']})", f"Yes\n(n={counts['Yes']})"])
# Save the plot
plt.xlabel('Reported PPI')
plt.ylabel('Cooperativity Index')
plt.tight_layout()
plt.savefig('ppi_reported_ppi.pdf')
plt.close()



#---------------------------------
# 2. Plot distribution of CI split by:
# 'Mediator','TFIIB','TFIID','TFIIE','TFIIF', 'TFIIH','SWI_SNF', 'POLII'
#---------------------------------

for col in ['Mediator','TFIIB','TFIID','TFIIE','TFIIF', 'TFIIH','SWI_SNF', 'POLII']:
    # Perform pairwise Mann-Whitney U tests
    categories = df[col].unique()
    p_values = []
    for cat1, cat2 in combinations(categories, 2):
        group1 = df[df[col] == cat1]["cooperativity_index"]
        group2 = df[df[col] == cat2]["cooperativity_index"]
        stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
        p_values.append({"cat1": cat1, "cat2": cat2, "p_value": p_value})
    # Convert to DataFrame for easier handling
    p_values_df = pd.DataFrame(p_values)
    # Count the number of observations in each category
    category_counts = df[col].value_counts()
    # Plot violin plot
    plt.figure(figsize=(6, 6))
    sns.violinplot(data=df, x=col, y="cooperativity_index", cut=0)
    # Annotate plot with p-values
    y_max = df["cooperativity_index"].max() + 0.1  # Adjust for annotation placement
    annotation_offset = 0.1  # Increase spacing between annotation lines
    line_offset = 0.02       # Adjust vertical height of lines
    for i, row in p_values_df.iterrows():
        cat1, cat2, p_value = row["cat1"], row["cat2"], row["p_value"]
        x1, x2 = categories.tolist().index(cat1), categories.tolist().index(cat2)
        y = y_max + i * annotation_offset  # Increase vertical spacing
        plt.plot([x1, x1, x2, x2], [y, y + line_offset, y + line_offset, y], lw=1.5, color="black")  # Draw annotation lines
        plt.text((x1 + x2) / 2, y + line_offset + 0.01, f"p={p_value:.3e}", ha="center")  # Add p-value
    # Add counts to x-axis labels
    xticks_labels = [f"{cat}\n(n={category_counts[cat]})" for cat in categories]
    plt.xticks(ticks=range(len(categories)), labels=xticks_labels)
    # Finalize and save the plot
    plt.xlabel(col)
    plt.ylabel("Cooperativity Index")
    plt.tight_layout()
    plt.savefig(f"ppi_{col}.pdf")
    plt.close()








