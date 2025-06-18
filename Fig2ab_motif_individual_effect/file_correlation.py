import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42




TRACK = "cage"

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
enhancers_merged = pd.merge(enhancers_hepg2_data, enhancers_k562_data, on="protein", 
                            suffixes=("_hepg2", "_k562"), how="inner")
promoters_merged = pd.merge(promoters_hepg2_data, promoters_k562_data, on="protein", 
                            suffixes=("_hepg2", "_k562"), how="inner")

# Extract relevant dstat_isa_cage_activity columns
enhancers_hepg2_dstat = enhancers_merged[f"dstat_isa_{TRACK}_activity_hepg2"]
enhancers_k562_dstat = enhancers_merged[f"dstat_isa_{TRACK}_activity_k562"]
promoters_hepg2_dstat = promoters_merged[f"dstat_isa_{TRACK}_activity_hepg2"]
promoters_k562_dstat = promoters_merged[f"dstat_isa_{TRACK}_activity_k562"]

# Prepare data dictionary for dstat values only
data_dstat = {
    "Enhancer HepG2": enhancers_hepg2_dstat,
    "Enhancer K562": enhancers_k562_dstat,
    "Promoter HepG2": promoters_hepg2_dstat,
    "Promoter K562": promoters_k562_dstat,
}

# Define the four comparisons as tuples: (x_category, y_category, plot_title)
comparisons = [
    ("Promoter HepG2","Enhancer HepG2", "HepG2: promoter vs enhancer"),
    ("Promoter K562", "Enhancer K562", "K562: promoter vs enhancer"),
    ("Promoter HepG2", "Promoter K562", "Promoter: HepG2 vs K562"),
    ("Enhancer HepG2", "Enhancer K562", "Enhancer: HepG2 vs K562"),
]

# Create subplot layout

fig, axes = plt.subplots(2, 2, figsize=(100/25.4, 100/25.4), dpi=300)
axes = axes.flatten()

for ax, (x_cat, y_cat, title) in zip(axes, comparisons):
    # Scatter plot using the dstat columns
    sns.scatterplot(x=data_dstat[x_cat], y=data_dstat[y_cat], ax=ax, alpha=0.6, s=5, color='black')
    # Compute the Pearson correlation coefficient and annotate the plot
    corr = data_dstat[x_cat].corr(data_dstat[y_cat])
    ax.annotate(f"r = {corr:.2f}", xy=(0.5, 0.9), xycoords='axes fraction',
                ha='center', fontsize=7, color='black')
    
    # Set titles and axis labels
    ax.set_title(title, fontsize=7)
    ax.set_xlabel(x_cat, fontsize=7)
    ax.set_ylabel(y_cat, fontsize=7)
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.6, 0.6)
    ax.tick_params(axis='both', which='major', labelsize=7)
    
# Adjust layout and add overall title if desired
plt.tight_layout()
plt.savefig(f"file_dstat_comparisons_{TRACK}.pdf")
plt.close()
