# analyze_dhs.py

import pandas as pd
import numpy as np
from loguru import logger
import sys
import matplotlib.pyplot as plt

# Ensure custom module paths
sys.path.insert(1, "/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity
from plotting_utils import plot_combined
from plotting import format_text

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# ----------------------------------------------------
# Constants
# ----------------------------------------------------
prefix = "/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs"

REGION_ORDER = ["distal_ts", "distal_ti", "proximal_ts", "proximal_ti"]
NEIGHBOR_PAIRS = [
    ("distal_ts", "distal_ti"),
    ("distal_ti", "proximal_ts"),
    ("proximal_ts", "proximal_ti"),
]

COLOR_MAP = {
    "redundant_tf_count":    "#1f77b4",
    "synergistic_tf_count":  "#d62728",
    "intermediate_tf_count": "grey",
    "independent_tf_count":  "black"
}

CELL_LABEL = {
    "hepg2": "HepG2",
    "k562":  "K562"
}


# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------
def read_df_add_tf_coop(prefix, file_name):
    df = pd.read_csv(f"{prefix}_{file_name}.csv", index_col=0)

    if "hepg2" in file_name:
        cell_line   = "hepg2"
        thresh_redun = 0.48
        thresh_codep = 0.78
    elif "k562" in file_name:
        cell_line   = "k562"
        thresh_redun = 0.43
        thresh_codep = 0.80
    else:
        raise ValueError("Unknown cell line in file_name")

    df_coop = pd.read_csv(
        f"/isdata/alab/people/pcr980/DeepCompare/"
        f"Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_dhs.csv"
    )
    df_coop = assign_cooperativity(df_coop, 5, 0.95, thresh_redun, thresh_codep)
    df_coop.rename(columns={"protein2": "protein"}, inplace=True)

    df = pd.merge(df, df_coop, on="protein", how="inner")
    return df


def count_tf(df, file_name):
    counts = (
        df.groupby(["region", "cooperativity"])
          .size()
          .unstack(fill_value=0)
          .reset_index()
    )
    counts.rename(columns={
        "Synergistic":  "synergistic_tf_count",
        "Redundant":    "redundant_tf_count",
        "Intermediate": "intermediate_tf_count",
        "Independent":  "independent_tf_count"
    }, inplace=True)
    counts["dataset"] = file_name
    return counts


def assign_region_type(df):
    df["re"]        = df["dataset"].apply(lambda x: x.split("_")[0])
    df["cell_line"] = df["dataset"].apply(lambda x: x.split("_")[1])
    conditions = [
        (df["re"] == "proximal") & (df["tissue_invariance"] == "yes"),
        (df["re"] == "proximal") & (df["tissue_invariance"] == "no"),
        (df["re"] == "distal")   & (df["tissue_invariance"] == "yes"),
        (df["re"] == "distal")   & (df["tissue_invariance"] == "no"),
    ]
    choices = ["proximal_ti", "proximal_ts", "distal_ti", "distal_ts"]
    df["region_type"] = np.select(conditions, choices, default="distal_ts")
    df["region_type"] = pd.Categorical(df["region_type"], categories=REGION_ORDER, ordered=True)
    return df


def add_seq_info(df, file_name):
    df_seq_info = pd.read_csv(
        f"/isdata/alab/people/pcr980/DeepCompare/"
        f"Pd4_promoters_enhancers_and_featimp/dhs_{file_name}.tsv",
        sep=" "
    )
    df_seq_info["region"] = (
        df_seq_info["seqnames"] + ":" +
        df_seq_info["start"].astype(str) + "-" +
        df_seq_info["end"].astype(str)
    )
    return pd.merge(df, df_seq_info, on="region", how="inner")


# ----------------------------------------------------
# Main
# ----------------------------------------------------
if __name__ == "__main__":

    # ---- A. TF counts (per region) ----
    def preprocess_counts(prefix, file_name):
        df_raw = read_df_add_tf_coop(prefix, file_name)
        df_cnt = count_tf(df_raw, file_name)
        df_cnt = add_seq_info(df_cnt, file_name)
        df_cnt["dataset"] = file_name
        return df_cnt

    df_counts_all = pd.concat([
        preprocess_counts(prefix, f"{re}_{cl}")
        for re in ["proximal", "distal"]
        for cl in ["hepg2", "k562"]
    ])
    df_counts_all = assign_region_type(df_counts_all)
    df_counts_all["region_type"] = df_counts_all["region_type"].cat.rename_categories(format_text)

    # ---- B. Cooperativity index ----
    def preprocess_ci(prefix, file_name):
        df = read_df_add_tf_coop(prefix, file_name)
        df["dataset"] = file_name
        df = add_seq_info(df, file_name)
        df = df[~df["cooperativity_index"].isna()].reset_index(drop=True)
        return assign_region_type(df)

    df_ci_all = pd.concat([
        preprocess_ci(prefix, f"{re}_{cl}")
        for re in ["proximal", "distal"]
        for cl in ["hepg2", "k562"]
    ])
    df_ci_all["region_type"] = df_ci_all["region_type"].cat.rename_categories(format_text)

    # ---- C. Build combined 2x3 panel (per cell line) ----
    for cell in ["hepg2", "k562"]:
        df_cnt = df_counts_all[df_counts_all["cell_line"] == cell].reset_index(drop=True)
        df_ci  = df_ci_all[df_ci_all["cell_line"] == cell].reset_index(drop=True)
        df_int = df_ci[df_ci["cooperativity"] == "Intermediate"]

        # Compute shared xlim for synergy score panels
        xmin = min(df_ci["cooperativity_index"].min(), df_int["cooperativity_index"].min())
        xmax = max(df_ci["cooperativity_index"].max(), df_int["cooperativity_index"].max())
        shared_xlim = (xmin, xmax)

        panel_specs = [
            {"df": df_ci,  "x_col": "cooperativity_index",   "title": "All TFs",           "color": None,  "showfliers": True,  "box_alpha": 0.3, "xlim": shared_xlim},
            {"df": df_cnt, "x_col": "synergistic_tf_count",  "title": "Synergistic",       "color": COLOR_MAP["synergistic_tf_count"],  "showfliers": False, "box_alpha": 0.3},
            {"df": df_cnt, "x_col": "redundant_tf_count",    "title": "Redundant",         "color": COLOR_MAP["redundant_tf_count"],    "showfliers": False, "box_alpha": 0.3},
            {"df": df_int, "x_col": "cooperativity_index",   "title": "Intermediate TFs",  "color": None,  "showfliers": True,  "box_alpha": 0.3, "xlim": shared_xlim},
            {"df": df_cnt, "x_col": "intermediate_tf_count", "title": "Intermediate",      "color": COLOR_MAP["intermediate_tf_count"], "showfliers": False, "box_alpha": 0.3},
            {"df": df_cnt, "x_col": "independent_tf_count",  "title": "Independent",       "color": COLOR_MAP["independent_tf_count"],  "showfliers": False, "box_alpha": 0.3},
        ]

        plot_combined(
            panel_specs    = panel_specs,
            neighbor_pairs = NEIGHBOR_PAIRS,
            outfile        = f"dhs_{cell}.pdf",
            figsize        = (5.0, 2.5),  
            nrows          = 2,
            ncols          = 3,
            sharey         = True,
            wspace         = 0.7,
            hspace         = 0.6
        )




# nohup python3 analyze_dhs.py > analyze_dhs.log 2>&1 &















