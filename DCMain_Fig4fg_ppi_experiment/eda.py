import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np

from utils_ppi import read_pooled_found_tf_pairs

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity


from adjustText import adjust_text

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42



suffix = "pe"
df_bait_pooled = read_pooled_found_tf_pairs()

for bait in ["BACH1", "RFX5", "IKZF1", "MAFG", "RREB1"]:
    for mode in ["linearity_index", "cooperativity_index"]:
        df_coop = pd.read_csv(
            f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_{suffix}.csv"
        )
        df_coop = assign_cooperativity(df_coop, 1, 0.9, 0.3, 0.7)

        df_ppi = pd.read_csv(
            f"/isdata/alab/people/pcr980/DeepCompare/DCMain_Fig4abcd_tf_pair_cooperativity/Pd1_Petra_data/2025-04-29_s10_PublishedPPIandProtComplexes_k562_pe.txt",
            sep="\t",
        )
        df_ppi.rename(columns={"SWI/SNF": "SWI_SNF"}, inplace=True)
        df_coop = pd.merge(df_coop, df_ppi, on=["protein1", "protein2"], how="left")

        df_coop["Reported_PPI"] = df_coop["Reported_PPI"].fillna(0)
        df_coop["Reported_PPI"] = df_coop["Reported_PPI"].map({"No": 0, "Yes": 1})
        df_coop = df_coop.dropna(subset=["Reported_PPI"]).reset_index(drop=True)

        df_coop = df_coop[df_coop["protein2"] == bait].reset_index(drop=True)

        df_bait = df_bait_pooled[df_bait_pooled["protein2"] == bait].reset_index(drop=True)
        df_coop["found"] = df_coop["protein1"].isin(df_bait["protein1"])

        # Create group labels for y-axis
        def get_group(row):
            if row["Reported_PPI"] == 1 and row["found"]:
                return "Reported, found"
            elif row["Reported_PPI"] == 1:
                return "Reported, not found"
            elif row["found"]:
                return "Not Reported, found"
            else:
                return "Not Reported, not found"

        df_coop["group"] = df_coop.apply(get_group, axis=1)

        order = [
            "Reported, found",
            "Reported, not found",
            "Not Reported, found",
            "Not Reported, not found",
        ]

        plt.figure(figsize=(3, 1.9))
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)

        # Draw stripplot
        sns.stripplot(
            data=df_coop,
            x=mode,
            y="group",
            order=order,
            size=2.5,
            alpha=0.2,
            jitter=0.15,
            color="black",
            dodge=False,
        )

        # Draw median line for each group
        grouped = df_coop.groupby("group")[mode].median()
        for i, grp in enumerate(order):
            median = grouped.get(grp, None)
            if median is not None:
                plt.plot([median, median], [i - 0.2, i + 0.2], lw=2, color="red")

        # Annotate names for found proteins (optional)
        if mode == "cooperativity_index":
            np.random.seed(42)  # For reproducible jitter
            for _, row in df_coop[df_coop["found"]].iterrows():
                y_val = order.index(row["group"])
                jitter = np.random.uniform(-0.5, 0.5)
                plt.text(
                    row[mode],
                    y_val + jitter,
                    row["protein1"],
                    fontsize=5,
                    color="black",
                    ha="center"
                )

        dict_x_label = {
            "linearity_index": "Independence score",
            "cooperativity_index": "Synergy score"
        }

        plt.xlabel(dict_x_label[mode], fontsize=7)
        plt.ylabel(None)
        plt.title(bait, fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(ticks=range(4), labels=order, fontsize=5)
        plt.legend([], [], frameon=False)  # No legend
        plt.tight_layout()
        plt.savefig(f"eda_{mode}_{bait}_{suffix}.pdf")
        plt.close()
