import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['pdf.fonttype']=42



TRACK="sure"

# Load datasets
enhancers_hepg2_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_enhancers_hepg2.csv'
enhancers_k562_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_enhancers_k562.csv'
promoters_hepg2_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_promoters_hepg2.csv'
promoters_k562_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_promoters_k562.csv'

enhancers_hepg2_data = pd.read_csv(enhancers_hepg2_file)
enhancers_k562_data = pd.read_csv(enhancers_k562_file)
promoters_hepg2_data = pd.read_csv(promoters_hepg2_file)
promoters_k562_data = pd.read_csv(promoters_k562_file)

# Merge datasets by protein for inner join
enhancers_merged = pd.merge(enhancers_hepg2_data, enhancers_k562_data, on="protein", suffixes=("_hepg2", "_k562"), how="inner")
promoters_merged = pd.merge(promoters_hepg2_data, promoters_k562_data, on="protein", suffixes=("_hepg2", "_k562"), how="inner")

# Extract relevant columns for avg_isa_cage_activity
enhancers_hepg2 = enhancers_merged[f"avg_isa_{TRACK}_activity_hepg2"]
enhancers_k562 = enhancers_merged[f"avg_isa_{TRACK}_activity_k562"]
promoters_hepg2 = promoters_merged[f"avg_isa_{TRACK}_activity_hepg2"]
promoters_k562 = promoters_merged[f"avg_isa_{TRACK}_activity_k562"]

# Extract relevant columns for dstat_isa_cage_activity
enhancers_hepg2_dstat = enhancers_merged[f"dstat_isa_{TRACK}_activity_hepg2"]
enhancers_k562_dstat = enhancers_merged[f"dstat_isa_{TRACK}_activity_k562"]
promoters_hepg2_dstat = promoters_merged[f"dstat_isa_{TRACK}_activity_hepg2"]
promoters_k562_dstat = promoters_merged[f"dstat_isa_{TRACK}_activity_k562"]

# Prepare data dictionaries
data = {
    "Enhancer HepG2": enhancers_hepg2,
    "Enhancer K562": enhancers_k562,
    "Promoter HepG2": promoters_hepg2,
    "Promoter K562": promoters_k562,
}

data_dstat = {
    "Enhancer HepG2": enhancers_hepg2_dstat,
    "Enhancer K562": enhancers_k562_dstat,
    "Promoter HepG2": promoters_hepg2_dstat,
    "Promoter K562": promoters_k562_dstat,
}

# Create the 4x4 super plot
fig, axes = plt.subplots(4, 4, figsize=(20, 20))
fig.suptitle("4x4 Super Plot: Lower Diagonal (avg_isa_cage), Upper Diagonal (dstat_cage)", fontsize=20)

categories = list(data.keys())

for i, row_category in enumerate(categories):
    for j, col_category in enumerate(categories):
        ax = axes[i, j]

        if i == j:
            # Remove diagonals by leaving them blank
            ax.axis('off')
        elif i > j:
            # Plot lower diagonal (avg_isa_cage comparisons)
            sns.scatterplot(x=data[row_category], y=data[col_category], ax=ax, alpha=0.6, edgecolor='none')
            corr = data[row_category].corr(data[col_category])
            ax.annotate(f"r={corr:.2f}", xy=(0.5, 0.9), xycoords='axes fraction',
                        ha='center', fontsize=10, color='blue')
            ax.grid(False)  # Remove the grid
        else:
            # Plot upper diagonal (dstat_cage comparisons)
            sns.scatterplot(x=data_dstat[row_category], y=data_dstat[col_category], ax=ax, alpha=0.6, edgecolor='none')
            corr = data_dstat[row_category].corr(data_dstat[col_category])
            ax.annotate(f"r={corr:.2f}", xy=(0.5, 0.9), xycoords='axes fraction',
                        ha='center', fontsize=10, color='red')
            ax.grid(False)  # Remove the grid

        # Label rows and columns
        if j == 0:
            ax.set_ylabel(row_category, fontsize=12)
        else:
            ax.set_ylabel("")

        if i == 3:
            ax.set_xlabel(col_category, fontsize=12)
        else:
            ax.set_xlabel("")

# Adjust layout
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(f"file_correlation_{TRACK}.pdf")
plt.close()
