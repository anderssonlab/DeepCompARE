import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from scipy.stats import mannwhitneyu
from loguru import logger

for cell_line in ["hepg2", "k562"]:
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    axes = axes.flatten()
    
    plot_idx = 0
    for tpm_thresh in [0.1, 0.5, 1.0]:
        for nonlinearity_thresh in [0.01, 0.05, 0.1]:
            ax = axes[plot_idx]
            plot_idx += 1

            # Logging
            logger.info(f"Processing {cell_line} {tpm_thresh} {nonlinearity_thresh}")
            
            # Load data
            df_tf = pd.read_csv(f"Pd1_cis/tf_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.csv")
            df_tf = df_tf[df_tf["c_sum"] > 5].reset_index(drop=True)

            # Sort by cooperativity index
            df_tf["rank"] = df_tf["cooperativity_index"].rank(ascending=True)

            # Load factor lists
            universal_stripe_factors = pd.read_csv("/isdata/alab/people/pcr980/Resource/universal_stripe_factors.txt", sep='\t').iloc[:, 0].tolist()
            pioneer_factors = pd.read_csv("/isdata/alab/people/pcr980/Resource/pioneer_factor_list.txt", sep='\t').iloc[:, 0].tolist()

            df_usf = df_tf[df_tf["protein2"].isin(universal_stripe_factors)]
            df_pf = df_tf[df_tf["protein2"].isin(pioneer_factors)]

            # Scatter plot
            sns.scatterplot(x="rank", y="cooperativity_index", data=df_tf, s=5, color='black', ax=ax)
            ax.set_xlim(-40, df_tf.shape[0] + 40)
            ax.set_ylim(-0.1, 1.1)
            ax.set_xlabel("Rank")
            ax.set_ylabel("TF cooperativity index")

            # Add text labels
            texts = []
            for i in range(df_usf.shape[0]):
                texts.append(ax.text(df_usf.iloc[i]["rank"], df_usf.iloc[i]["cooperativity_index"], df_usf.iloc[i]["protein2"], color='#4169E1'))

            for i in range(df_pf.shape[0]):
                texts.append(ax.text(df_pf.iloc[i]["rank"], df_pf.iloc[i]["cooperativity_index"], df_pf.iloc[i]["protein2"], color='darkorange'))

            adjust_text(texts, ax=ax)

            # Add legend for text color
            ax.scatter([], [], color='#4169E1', label='Universal stripe factors')
            ax.scatter([], [], color='darkorange', label='Pioneer factors')
            ax.legend()

            # Mann-Whitney U test
            df_non_usf = df_tf[~df_tf["protein2"].isin(universal_stripe_factors)]
            _, p_usf = mannwhitneyu(df_usf["cooperativity_index"], df_non_usf["cooperativity_index"], alternative='less')

            df_non_pf = df_tf[~df_tf["protein2"].isin(pioneer_factors)]
            _, p_pf = mannwhitneyu(df_pf["cooperativity_index"], df_non_pf["cooperativity_index"], alternative='greater')

            ax.set_title(f"tpm={tpm_thresh}, nonlinearity={nonlinearity_thresh}\nUSF p={p_usf:.2e}, PF p={p_pf:.2e}")

    plt.tight_layout()
    plt.savefig(f"Plots/usf_pf_{cell_line}.png")
    plt.close()
