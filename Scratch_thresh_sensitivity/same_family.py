import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load family data
df_family = pd.read_csv("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2024_CORE_tf_family.csv")
df_family["ID"] = df_family["ID"].str.upper()

cell_lines = ["hepg2", "k562"]
tpm_thresholds = [0.1, 0.5, 1.0]
nonlinearity_thresholds = [0.01, 0.05, 0.1]

for cell_line in cell_lines:
    fig, axes = plt.subplots(len(tpm_thresholds), len(nonlinearity_thresholds), figsize=(15,15), sharex=True, sharey=True)
    fig.suptitle(f"Aggregated Violin Plots for {cell_line}", fontsize=16)

    for i, tpm_thresh in enumerate(tpm_thresholds):
        for j, nonlinearity_thresh in enumerate(nonlinearity_thresholds):
            ax = axes[i, j]

            # Load data
            df = pd.read_csv(f"Pd1_cis/tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.csv")
            df = df[df["c_sum"] > 1].reset_index(drop=True)

            # Merge with family data
            df = pd.merge(df, df_family, left_on="protein1", right_on="ID", how="inner")
            df.drop(columns=['AC', 'ID'], inplace=True)
            df.rename(columns={"tf_family": "family_1"}, inplace=True)
            df = pd.merge(df, df_family, left_on="protein2", right_on="ID", how="inner")
            df.drop(columns=['AC', 'ID'], inplace=True)
            df.rename(columns={"tf_family": "family_2"}, inplace=True)

            # Determine if same family
            df["same_family"] = df["family_1"] == df["family_2"]

            # Filter families with at least 20 pairs
            df_filtered = pd.DataFrame()
            for family, df_sub in df.groupby("family_2"):
                if df_sub["same_family"].sum() >= 20:
                    df_filtered = pd.concat([df_filtered, df_sub], axis=0)

            if not df_filtered.empty:
                # Calculate family sizes
                family_sizes = df_filtered["family_2"].value_counts().reset_index()
                family_sizes.columns = ["family_2", "family_size"]
                df_filtered = pd.merge(df_filtered, family_sizes, on="family_2", how="left")

                # Plot violin plot
                sns.violinplot(
                    data=df_filtered,
                    x="family_2",
                    y="cooperativity_index",
                    hue="same_family",
                    split=True,
                    cut=0,
                    ax=ax,
                    inner=None
                )
                ax.set_xticks(range(len(family_sizes)))
                ax.set_xticklabels(
                    [f"{family}\n(n={size})" for family, size in zip(family_sizes["family_2"], family_sizes["family_size"])],
                    rotation=90
                )
            else:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=12)
                ax.set_xticks([])
                ax.set_yticks([])

            ax.set_title(f"TPM={tpm_thresh}, NL={nonlinearity_thresh}")

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"same_family_{cell_line}.png")
    plt.close()
