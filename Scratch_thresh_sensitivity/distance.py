import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,mannwhitneyu
from adjustText import adjust_text

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label
from plotting import plot_violin_with_statistics



import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Iterate through each cell line
for cell_line in ["hepg2", "k562"]:
    # Create a PDF to save all plots in a single file
    with PdfPages(f'distance_vs_cooperativity_index_{cell_line}.pdf') as pdf:
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
        fig.suptitle(f'Distance vs Cooperativity Index Bin for {cell_line}', fontsize=16)
        
        # Iterate through tpm_thresh and nonlinearity_thresh combinations
        for i, tpm_thresh in enumerate([0.1, 0.5, 1.0]):
            for j, nonlinearity_thresh in enumerate([0.01, 0.05, 0.1]):
                # Read and process data
                df = pd.read_csv(f"Pd1_cis/tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.csv")
                df = df[df["c_sum"] > 1].reset_index(drop=True)
                df = bin_and_label(df, "cooperativity_index", [0, 0.3, 0.7, 1])
                df["cooperativity"] = np.where(
                    df["cooperativity_index"] > 0.7, 
                    "codependent", 
                    np.where(df["cooperativity_index"] < 0.3, "redundant", "")
                )
                
                # Save individual plot to a temporary file
                temp_output_file = f"temp_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.png"
                plot_violin_with_statistics(
                    df=df,
                    x_col="cooperativity_index_bin",
                    y_col="distance",
                    x_label="Cooperativity Index Bin",
                    y_label="Median Distance (bp)",
                    title=f"TPM={tpm_thresh}, Nonlin={nonlinearity_thresh}",
                    output_file=temp_output_file  # Save the individual plot
                )
                
                # Load the saved plot and insert it into the grid
                img = plt.imread(temp_output_file)
                axes[i, j].imshow(img)
                axes[i, j].axis('off')  # Hide axes since we're displaying an image
            
        # Save the combined figure to the PDF
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)
