import pandas as pd
import numpy as np
from loguru import logger
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from stat_tests import bin_and_label, calc_or_by_various_thresholds
from plotting import format_text




import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


def plot_or(df, color_mapping, n_mapping, title, out_path):
    """
    df: output of odds_ratio_one_df, contains 'or', 'pval', 'ci_low', 'ci_high', 'threshold'
    """
    # Format the title for display
    title = format_text(title)
    plt.figure(figsize=(2.65, 2.2))
    # thinner frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    # determine transparency
    df["alphas"] = df["pval"].apply(lambda x: 1.0 if x < 0.05 else 0.1)
    # x is the order of bin in categorical
    df["bin"] = df.index
    df["x"] = df["bin"].cat.codes
    jitter = 0
    for threshold, df_subset in df.groupby(["threshold"]):
        plt.scatter(df_subset["x"] + jitter,
                    df_subset["or"],
                    color=color_mapping[threshold],
                    alpha=df_subset['alphas'],
                    s=7)
        for _, row in df_subset.iterrows():
            plt.errorbar(
                row["x"] + jitter,
                row["or"],
                yerr=[[row["or"] - row["ci_low"]], [row["ci_high"] - row["or"]]],
                fmt='none',
                color=color_mapping[threshold],
                capsize=0,
                alpha=row["alphas"],
                markeredgewidth=0.5,
                elinewidth=1)
        plt.scatter(
            [], [],  # Invisible data
            color=color_mapping[threshold],
            label=f"AF>{threshold} (n={n_mapping[threshold]})", s=10)
        jitter += 0.15
    plt.xlabel("Predicted loss-of-function effect size", fontsize=7)
    plt.ylabel("Odds ratio", fontsize=7)
    plt.title(title, fontsize=7)
    plt.axhline(y=1, color='black', linestyle=':', linewidth=0.5)
    plt.legend(fontsize=5)
    plt.xticks(df["x"], df["bin"], rotation=30, fontsize=5)
    plt.yticks(fontsize=5)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()




#----------------------
# Read data
#----------------------
threshold_list = [0.0001, 0.001, 0.01, 0.05]
colors_list = ['#808080', '#B5BD00', '#1E90FF', '#000080']
color_mapping = dict(zip(threshold_list, colors_list))

# Maps dataset names to file paths
dataset_to_file = {
    "promoters_hepg2": "Pd1_maf_with_effect_size/maf_with_effect_size_promoters_hepg2.csv",
    "promoters_k562": "Pd1_maf_with_effect_size/maf_with_effect_size_promoters_k562.csv",
    "enhancers_hepg2": "Pd1_maf_with_effect_size/maf_with_effect_size_enhancers_hepg2.csv",
    "enhancers_k562": "Pd1_maf_with_effect_size/maf_with_effect_size_enhancers_k562.csv",
}

track_num_2_name = {
    0: "CAGE",
    1: "CAGE",
    2: "DNase",
    3: "DNAse",
    4: "STARR",
    5: "STARR",
    6: "SuRE",
    7: "SuRE",
}

for dataset, file_path in dataset_to_file.items():
    df = pd.read_csv(file_path, header=None, index_col=0)
    df["dataset"] = dataset
    # Rename columns only once
    df.columns = ["chromosome", "start", "end", "ID", "REF", "ALT", "AF", "Name", "Score", "Strand"] + \
                 [f"track_{i}" for i in range(16)] + ["dataset"]
    for track_num in get_track_num(dataset):
        logger.info(f"Processing {dataset}, track {track_num}")
        # Filter and preprocess
        df_filtered = df[df[f"track_{track_num}"] < 0].copy()
        df_filtered[f"track_{track_num}"] = df_filtered[f"track_{track_num}"].abs()
        df_filtered = bin_and_label(df_filtered, f"track_{track_num}", [0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])
        # Calculate OR
        df_or, n_list = calc_or_by_various_thresholds(df_filtered, "AF", threshold_list, "larger", f"track_{track_num}_bin")
        # Save results and plot
        df_or.to_csv(f"Tables_or/or_{dataset}_track{track_num}.csv", index=False)
        n_mapping = dict(zip(threshold_list, n_list))
        plot_or(df_or, color_mapping, n_mapping, f"{dataset}, {track_num_2_name[track_num]} effect", f"or_{dataset}_track{track_num}.pdf")



# nohup python3 plot_or.py > plot_or.out &
