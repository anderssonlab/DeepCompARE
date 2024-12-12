import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

for cell_type in ["hepg2","k562"]:
    for tpm_thresh in [0.1,0.5,1.0]:
        for nonlinearity_thresh in [0.01,0.05,0.1]:
            # 1. histogram of tf_pair ci
            df=pd.read_csv(f"Pd1_cis/tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_type}.csv")
            df=df[df["c_sum"]>1].reset_index(drop=True)
            sns.histplot(df['cooperativity_index'], bins=100)
            plt.xlabel('TF pair cooperativity index')
            plt.ylabel('Frequency')
            plt.title(f'tpm_thresh={tpm_thresh}, nonlinearity_thresh={nonlinearity_thresh} ({cell_type})')
            plt.savefig(f'hist_tf_pair_ci_{tpm_thresh}_{nonlinearity_thresh}_{cell_type}.png')
            plt.close()

