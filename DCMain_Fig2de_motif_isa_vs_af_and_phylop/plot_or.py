import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from stat_tests import bin_and_label, calc_or_by_various_thresholds
from plotting import format_text

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
#-------------------
# Helper functions
#-------------------
def read_file(file_suffix):
    file_path = f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv"
    df = pd.read_csv(file_path, index_col=0)
    df["max_af"] = df["gnomad_af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["max_241way"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["dataset"] = file_suffix
    df["cell_type"] = file_suffix.split("_")[1]
    return df

def plot_or(df, color_mapping, n_mapping, feature_label, title, out_path):
    plt.figure(figsize=(2.65, 2.2))
    ax = plt.gca()
    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_linewidth(0.5)

    df["alphas"] = df["pval"].apply(lambda x: 1.0 if x < 0.05 else 0.1)
    df["bin"] = df.index
    df["x"] = df["bin"].cat.codes
    jitter = 0

    for threshold, df_subset in df.groupby("threshold"):
        ax.scatter(df_subset["x"] + jitter, df_subset["or"],
                   color=color_mapping[threshold], alpha=df_subset["alphas"], s=7)
        for _, row in df_subset.iterrows():
            ax.errorbar(row["x"] + jitter, row["or"],
                        yerr=[[row["or"] - row["ci_low"]], [row["ci_high"] - row["or"]]],
                        fmt='none', color=color_mapping[threshold], alpha=row["alphas"],
                        capsize=0, markeredgewidth=0.5, elinewidth=1)
        ax.scatter([], [], color=color_mapping[threshold],
                   label=f"max({feature_label})>{threshold} (n={n_mapping[threshold]})", s=10)
        jitter += 0.15

    ax.set_xlabel("Motif-level ISA", fontsize=7)
    ax.set_ylabel("Odds ratio", fontsize=7)
    ax.axhline(y=1, color='black', linestyle=':', linewidth=0.5)
    ax.legend(fontsize=5)
    ax.set_title(format_text(title), fontsize=7)
    ax.set_xticks(df["x"])
    ax.set_xticklabels(df["bin"], rotation=30, fontsize=5)
    ax.tick_params(axis='y', labelsize=5)

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

#-------------------
# Main logic
#-------------------

track_num_2_name = {
    0: "CAGE", 1: "CAGE", 2: "DNase", 3: "DNase",
    4: "STARR", 5: "STARR", 6: "SuRE", 7: "SuRE",
}

# Plot max_af
af_threshold_list = [0.0001, 0.001, 0.01, 0.05]
af_colors = ["#808080", "#FFC0CB", "#C7B6EA", "#9400D3"]

# Plot max_241way
phylop_threshold_list = [0, 1, 2, 3]
phylop_colors = ["#808080", "#FFD700", "#FFA500", "#D98880"]




for dataset in ["promoters_hepg2", "enhancers_hepg2", "promoters_k562", "enhancers_k562"]:
    df = read_file(dataset)
    for track_num in get_track_num(dataset):
        logger.info(f"Processing {dataset}, track {track_num}")
        track_col = f"isa_track{track_num}"
        binned_col = f"{track_col}_bin"
        df_filtered = df[df[track_col] > 0].reset_index(drop=True)
        df_binned = bin_and_label(df_filtered, track_col, [0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])
        # Plot max_af odds ratio
        df_or_af, n_list_af = calc_or_by_various_thresholds(df_binned, "max_af",
                                                            af_threshold_list,
                                                            "larger",
                                                            binned_col)
        df_or_af.to_csv(f"Tables_af/or_{dataset}_track{track_num}.csv")
        plot_or(
            df_or_af,
            dict(zip(af_threshold_list, af_colors)),
            dict(zip(af_threshold_list, n_list_af)),
            "AF",
            f"{format_text(dataset)}, {track_num_2_name[track_num]} effect",
            f"af_or_{dataset}_track{track_num}.pdf"
        )

        # Plot phylop odds ratio
        df_or_phylo, n_list_phylo = calc_or_by_various_thresholds(df_binned, "max_241way",
                                                                  phylop_threshold_list,
                                                                  "larger",
                                                                  binned_col)
        df_or_phylo.to_csv(f"Tables_phylop/or_{dataset}_track{track_num}.csv")
        plot_or(
            df_or_phylo,
            dict(zip(phylop_threshold_list, phylop_colors)),
            dict(zip(phylop_threshold_list, n_list_phylo)),
            "phylop",
            f"{format_text(dataset)}, {track_num_2_name[track_num]} effect",
            f"phylop_or_{dataset}_track{track_num}.pdf"
        )

# nohup python3 plot_or.py > plot_or.out &