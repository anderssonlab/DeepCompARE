# analyze_e1e2p1p2.py

import pandas as pd
import sys

sys.path.insert(1, "/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity
from plotting_utils import plot_combined
from plotting import format_text

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

prefix = "Pd2_E1E2P1P2_motif_info/motif_info_thresh_500"
REGION_ORDER = ["E1", "E2", "P2", "P1"]
NEIGHBOR_PAIRS = [("E1", "E2"), ("E2", "P2"), ("P2", "P1")]
COLOR_MAP = {
    "redundant_tf_count": "#1f77b4",
    "synergistic_tf_count": "#d62728",
    "intermediate_tf_count": "grey",
    "independent_tf_count": "black"
}

def read_df_add_tf_coop(prefix, file_name):
    df = pd.read_csv(f"{prefix}_{file_name}_k562.csv", index_col=0)
    df_coop = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_k562_dhs.csv")
    df_coop = assign_cooperativity(df_coop, 5, 0.95, 0.43, 0.80)
    df_coop.rename(columns={"protein2": "protein"}, inplace=True)
    df = pd.merge(df, df_coop, on="protein", how="inner")
    df["region_type"] = file_name
    return df

def count_tf(df, file_name):
    counts = (
        df.groupby(["region", "cooperativity"]).size()
        .unstack(fill_value=0).reset_index()
    )
    counts.rename(columns={
        "Synergistic": "synergistic_tf_count",
        "Redundant": "redundant_tf_count",
        "Intermediate": "intermediate_tf_count",
        "Independent": "independent_tf_count"
    }, inplace=True)
    counts["region_type"] = file_name
    return counts

def prepare_data(prefix):
    df_p1 = read_df_add_tf_coop(prefix, "P1")
    df_p2 = read_df_add_tf_coop(prefix, "P2")
    df_e1 = read_df_add_tf_coop(prefix, "E1")
    df_e2 = read_df_add_tf_coop(prefix, "E2")

    df_counts = pd.concat([
        count_tf(df_p1, "P1"), count_tf(df_p2, "P2"),
        count_tf(df_e1, "E1"), count_tf(df_e2, "E2")
    ], axis=0)
    df_counts["region_type"] = pd.Categorical(df_counts["region_type"], categories=REGION_ORDER, ordered=True)

    df_full = pd.concat([df_p1, df_p2, df_e1, df_e2], axis=0)
    df_full["region_type"] = pd.Categorical(df_full["region_type"], categories=REGION_ORDER, ordered=True)

    df_intermediate = df_full[df_full["cooperativity"] == "Intermediate"].reset_index(drop=True)
    return df_counts, df_full, df_intermediate

if __name__ == "__main__":
    df_counts, df_full, df_intermediate = prepare_data(prefix)

    for df in (df_counts, df_full, df_intermediate):
        df["region_type"] = df["region_type"].cat.rename_categories(lambda s: format_text(s))

    xmin = min(df_full["cooperativity_index"].min(), df_intermediate["cooperativity_index"].min())
    xmax = max(df_full["cooperativity_index"].max(), df_intermediate["cooperativity_index"].max())
    shared_xlim = (xmin, xmax)

    panel_specs = [
        {"df": df_full,         "x_col": "cooperativity_index",   "title": "All TFs",         "color": None, "showfliers": True, "box_alpha": 0.3, "xlim": shared_xlim},
        {"df": df_counts,       "x_col": "synergistic_tf_count",  "title": "Synergistic",     "color": COLOR_MAP["synergistic_tf_count"], "showfliers": True, "box_alpha": 0.3},
        {"df": df_counts,       "x_col": "redundant_tf_count",    "title": "Redundant",       "color": COLOR_MAP["redundant_tf_count"], "showfliers": True, "box_alpha": 0.3},
        {"df": df_intermediate, "x_col": "cooperativity_index",   "title": "Intermediate TFs", "color": None, "showfliers": True, "box_alpha": 0.3, "xlim": shared_xlim},
        {"df": df_counts,       "x_col": "intermediate_tf_count", "title": "Intermediate",    "color": COLOR_MAP["intermediate_tf_count"], "showfliers": True, "box_alpha": 0.3},
        {"df": df_counts,       "x_col": "independent_tf_count",  "title": "Independent",     "color": COLOR_MAP["independent_tf_count"], "showfliers": True, "box_alpha": 0.3},
    ]

    plot_combined(
        panel_specs=panel_specs,
        neighbor_pairs=NEIGHBOR_PAIRS,
        outfile="e1e2p1p2.pdf",
        figsize=(4.2, 2.5),
        nrows=2,
        ncols=3,
        sharey=True,
        wspace=0.7,
        hspace=0.6
    )
