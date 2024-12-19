import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# Iterate through each cell type
for cell_type in ["hepg2", "k562"]:
    fig, axes = plt.subplots(3, 3, figsize=(15, 15), sharex=True, sharey=True)
    fig.suptitle(f'Histograms of TF Pair Cooperativity Index for {cell_type}', fontsize=16)
    
    # Iterate through tpm_thresh and nonlinearity_thresh combinations
    for i, tpm_thresh in enumerate([0.1, 0.5, 1.0]):
        for j, nonlinearity_thresh in enumerate([0.01, 0.05, 0.1]):
            # Read data
            df = pd.read_csv(f"Pd1_cis/tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_type}.csv")
            df = df[df["c_sum"] > 1].reset_index(drop=True)
            
            # Create histogram in the appropriate subplot
            ax = axes[i, j]
            sns.histplot(df['cooperativity_index'], bins=100, ax=ax)
            ax.set_title(f'TPM={tpm_thresh}, Nonlin={nonlinearity_thresh}')
            ax.set_xlabel('Cooperativity Index')
            ax.set_ylabel('Frequency')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f'hist_ci_{cell_type}.png')
    plt.close()
