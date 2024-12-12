import pandas as pd


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics




#---------------------------------
# Read and merge data
#---------------------------------
for cell_line in ["hepg2", "k562"]:
    for tpm_thresh in [0.1,0.5,1.0]:
        for nonlinearity_thresh in [0.01,0.05,0.1]:
            df_coop = pd.read_csv(f"Pd1_cis/tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.csv")
            df_coop = df_coop[df_coop["c_sum"] > 1].reset_index(drop=True)
            df_ppi = pd.read_csv("/isdata/alab/people/pcr980/Resource/2024-10-24_s1_PublishedPPIannotation_withProtComplexes.txt", sep='\t')
            df_ppi.rename(columns={'SWI/SNF_lowerTh': 'SWI_SNF_lowerTh', 'SWI/SNF': 'SWI_SNF'}, inplace=True)
            df_ppi.drop(columns=['protPair', 'coopIndex_avg', 'sum_ci_avg'], inplace=True)
            df_ppi2 = df_ppi.copy()
            df_ppi2.rename(columns={'protein1': 'protein2', 'protein2': 'protein1'}, inplace=True)
            df_ppi = pd.concat([df_ppi, df_ppi2], axis=0).reset_index(drop=True)

            df = pd.merge(df_coop, df_ppi, left_on=['protein1', 'protein2'], right_on=['protein1', 'protein2'], how='inner')
            # Plot 1: Reported_PPI
            plot_violin_with_statistics(
                df=df,
                x_col="Reported_PPI",
                y_col="cooperativity_index",
                x_label="Reported PPI",
                y_label="Cooperativity Index",
                title=f"tpm_thresh={tpm_thresh}, nonlinearity_thresh={nonlinearity_thresh} ({cell_line})",
                output_file=f"ppi_reported_ppi_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.pdf"
            )

            # Plot 2: Protein complexes
            for col in ['Mediator', 'SWI_SNF']:
                plot_violin_with_statistics(
                    df=df,
                    x_col=col,
                    y_col="cooperativity_index",
                    x_label=col,
                    y_label="Cooperativity Index",
                    title=f"tpm_thresh={tpm_thresh}, nonlinearity_thresh={nonlinearity_thresh} ({cell_line})",
                    output_file=f"ppi_{col}_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.pdf"
                )