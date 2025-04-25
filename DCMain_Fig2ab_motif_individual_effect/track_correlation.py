import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


AGG_METHOD="dstat"

# Load the HepG2 dataset
hepg2_file_path = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_hepg2_pe.csv'
hepg2_data = pd.read_csv(hepg2_file_path, index_col=0)

# Load the K562 dataset
k562_file_path = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_k562_pe.csv'
k562_data = pd.read_csv(k562_file_path, index_col=0)

# Define the relevant columns
columns = [
    f"{AGG_METHOD}_isa_cage_activity",
    f"{AGG_METHOD}_isa_dhs_activity",
    f"{AGG_METHOD}_isa_starr_activity",
    f"{AGG_METHOD}_isa_sure_activity",
]

# Select subsets for HepG2 and K562
hepg2_subset = hepg2_data[columns]
k562_subset = k562_data[columns]

# merge on index
merged_subset=pd.merge(hepg2_subset, k562_subset, left_index=True, right_index=True, suffixes=('_hepg2', '_k562'))

# Create the scatter plot matrix
# Create the scatter plot matrix with requested updates
fig, axes = plt.subplots(4, 4, figsize=(170/25.4, 170/25.4), dpi=300, sharex='col', sharey='row')
# fig.suptitle("Track correlation: HepG2 (Lower) and K562 (Upper)", fontsize=7)

# Iterate through each pair of columns to plot
for i, col1 in enumerate(columns):
    for j, col2 in enumerate(columns):
        ax = axes[i, j]
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.grid(False)  # Remove the grid
        if i == j:
            # Plot HepG2 vs K562 on the diagonal, remove grid, and add Pearson r
            sns.scatterplot(x=hepg2_subset[col1], y=k562_subset[col1], ax=ax, alpha=0.6, edgecolor='none', color='#33b233',s=5)
            corr = hepg2_subset[col1].corr(k562_subset[col1])  # Pearson correlation for HepG2 vs K562
            ax.annotate(f"r={corr:.2f}", xy=(0.5, 0.9), xycoords='axes fraction',
                        ha='center', fontsize=5, color='#33b233')
        elif i > j:
            # Plot HepG2 on the lower diagonal
            sns.scatterplot(x=hepg2_subset[col1], y=hepg2_subset[col2], ax=ax, alpha=0.6, edgecolor='none', color='#d93333',s=5)
            # Calculate Pearson correlation for HepG2
            corr = hepg2_subset[col1].corr(hepg2_subset[col2])
            ax.annotate(f"r={corr:.2f}", xy=(0.5, 0.9), xycoords='axes fraction',
                        ha='center', fontsize=5, color='#d93333')
        else:
            # Plot K562 on the upper diagonal
            sns.scatterplot(x=k562_subset[col1], y=k562_subset[col2], ax=ax, alpha=0.6, edgecolor='none', color='#3366cc',s=5)
            # Calculate Pearson correlation for K562
            corr = k562_subset[col1].corr(k562_subset[col2])
            ax.annotate(f"r={corr:.2f}", xy=(0.5, 0.9), xycoords='axes fraction',
                        ha='center', fontsize=5, color='#3366cc')
        #
        # Retain y-labels for the first column
        if j == 0:
            label=col1.split("_")[2].upper()
            ax.set_ylabel(label, fontsize=7)
        else:
            ax.set_ylabel("")
        #
        # Retain x-labels for the last row
        if i == 3:
            label=col2.split("_")[2].upper()
            ax.set_xlabel(label, fontsize=7)
        else:
            ax.set_xlabel("")


# Adjust spacing
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(f'track_correlation_{AGG_METHOD}.pdf')
plt.close()

