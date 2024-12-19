import pandas as pd
from matplotlib import pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics


import pandas as pd
import matplotlib.pyplot as plt


# Load PPI data
df_ppi = pd.read_csv("/isdata/alab/people/pcr980/Resource/2024-10-24_s1_PublishedPPIannotation_withProtComplexes.txt", sep='\t')
df_ppi.rename(columns={'SWI/SNF_lowerTh': 'SWI_SNF_lowerTh', 'SWI/SNF': 'SWI_SNF'}, inplace=True)
df_ppi.drop(columns=['protPair', 'coopIndex_avg', 'sum_ci_avg'], inplace=True)
df_ppi2 = df_ppi.copy()
df_ppi2.rename(columns={'protein1': 'protein2', 'protein2': 'protein1'}, inplace=True)
df_ppi = pd.concat([df_ppi, df_ppi2], axis=0).reset_index(drop=True)

# Aggregated plot function
def aggregate_violin_plots_with_statistics(cell_line, tpm_thresholds, nonlinearity_thresholds):
    fig, axes = plt.subplots(len(tpm_thresholds), len(nonlinearity_thresholds), figsize=(15, 15))
    axes = axes.flatten()

    for idx, (tpm_thresh, nonlinearity_thresh) in enumerate(
            [(tpm, nl) for tpm in tpm_thresholds for nl in nonlinearity_thresholds]):
        ax = axes[idx]

        # Load cooperativity data
        df_coop = pd.read_csv(f"Pd1_cis/tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.csv")
        df_coop = df_coop[df_coop["c_sum"] > 1].reset_index(drop=True)

        # Merge with PPI data
        df = pd.merge(df_coop, df_ppi, left_on=['protein1', 'protein2'], right_on=['protein1', 'protein2'], how='inner')

        # Skip empty data
        if df.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        # Generate a temporary output file name for the subplot
        temp_output_file = f"temp_plot_{cell_line}_tpm{tpm_thresh}_nl{nonlinearity_thresh}.png"

        # Use your plot_violin_with_statistics function
        plot_violin_with_statistics(
            df=df,
            x_col="Reported_PPI",
            y_col="cooperativity_index",
            x_label="Reported PPI",
            y_label="Cooperativity Index",
            title=f"TPM={tpm_thresh}, NL={nonlinearity_thresh}",
            output_file=temp_output_file
        )

        # Load the saved plot and add it to the subplot
        img = plt.imread(temp_output_file)
        ax.imshow(img)
        ax.axis('off')

    # Final adjustments and save the aggregated figure
    plt.tight_layout()
    plt.suptitle(f"Aggregated Violin Plots for {cell_line}", fontsize=16, y=1.02)
    aggregated_output_file = f"ppi_{cell_line}.png"
    plt.savefig(aggregated_output_file)
    plt.close()

# Generate aggregated plots for each cell line
for cell_line in ["hepg2", "k562"]:
    aggregate_violin_plots_with_statistics(cell_line, [0.1, 0.5, 1.0], [0.01, 0.05, 0.1])
