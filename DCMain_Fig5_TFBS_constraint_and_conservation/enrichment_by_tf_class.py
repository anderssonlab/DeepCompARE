import pandas as pd
import numpy as np
from pybedtools import BedTool
from matplotlib import pyplot as plt
import seaborn as sns
import warnings


import sys
sys.path.append("/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import fisher_exact_with_ci






def add_tf_codependency(df,suffix):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{suffix}.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{suffix}.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df


def count_overlaps(feature, regions):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")  
        print(feature.intersect(regions).to_dataframe())
        return feature.intersect(regions).to_dataframe().shape[0]





cell_line="k562"
dict_track_num={"hepg2":0, "k562":1}

# Read BED files into BedTool objects
cds=BedTool("cds.bed")
trait = BedTool("/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions/whole_blood_eQTLs_PIP_01_non_coding.bed")
snp = BedTool("/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions/background_SNPs_non_coding.bed")



dhs_proximal=pd.read_csv(f"Pd1_dnase_annotated/dhs_proximal_{cell_line}.csv")
dhs_distal=pd.read_csv(f"Pd1_dnase_annotated/dhs_distal_{cell_line}.csv")
dhs=pd.concat([dhs_proximal, dhs_distal], axis=0, ignore_index=True)
dhs = BedTool.from_dataframe(dhs)


# Load TFBS DHS data
df_proximal = pd.read_csv(f"motif_info_thresh_500_dhs_proximal_{cell_line}.csv")
# df_distal = pd.read_csv(f"motif_info_thresh_500_dhs_distal_{cell_line}.csv")


df_distal_1=pd.read_csv("motif_info_thresh_500_dhs_distal_k562_v1.csv",nrows=938528)
df_distal_2=pd.read_csv("motif_info_thresh_500_dhs_distal_k562_v1.csv",skiprows=938529,header=None)
df_distal_2.drop(df_distal_2.columns[29], axis=1, inplace=True)
df_distal_2.columns=df_distal_1.columns
df_distal=pd.concat([df_distal_1,df_distal_2], axis=0, ignore_index=True)



df_proximal = add_tf_codependency(df_proximal,cell_line)
df_distal = add_tf_codependency(df_distal,cell_line)
# remove 0th column
df_proximal = df_proximal.iloc[:,1:]
df_distal = df_distal.iloc[:,1:]
# convert to BedTool
gr_proximal = BedTool.from_dataframe(df_proximal)
gr_distal = BedTool.from_dataframe(df_distal)


# Remove overlaps with CDS
gr_proximal = gr_proximal.intersect(cds, v=True)
gr_distal = gr_distal.intersect(cds, v=True)
dhs = dhs.intersect(cds, v=True)

# Calculate overlaps as base counts

n_total_background = count_overlaps(snp, dhs)
n_total_trait = count_overlaps(trait, dhs)



def calc_or_by_threshold(df_tfbs, bed_trait, snp, n_total_trait, n_total_background):
    df_res = []
    for threshold in np.arange(0, 0.5, 0.02):
        df_tfbs_above_thresh = df_tfbs[df_tfbs[f"isa_track{dict_track_num[cell_line]}"] > threshold].reset_index(drop=True)
        bed_tfbs_above_thresh = BedTool.from_dataframe(df_tfbs_above_thresh)
        n_tfbs_trait = count_overlaps(bed_tfbs_above_thresh, bed_trait)
        n_tfbs_background= count_overlaps(bed_tfbs_above_thresh, snp)
        odds_ratio, pval, ci_low, ci_high = fisher_exact_with_ci(np.array([[n_tfbs_trait, n_total_trait - n_tfbs_trait], [n_tfbs_background, n_total_background - n_tfbs_background]]))
        df_res.append({"threshold": threshold, "odds_ratio": odds_ratio, "pvalue": pval, "ci_low": ci_low, "ci_high": ci_high})
    return pd.DataFrame(df_res)




df_proximal = gr_proximal.to_dataframe(names=df_proximal.columns)
df_tfbs_proximal_redundant = df_proximal[df_proximal["cooperativity"]=="redundant"].reset_index(drop=True)
df_tfbs_proximal_codependent = df_proximal[df_proximal["cooperativity"]=="codependent"].reset_index(drop=True)

df_res_proximal_redundant = calc_or_by_threshold(df_tfbs_proximal_redundant, trait, snp, n_total_trait, n_total_background)
df_res_proximal_codependent = calc_or_by_threshold(df_tfbs_proximal_codependent, trait, snp, n_total_trait, n_total_background)



df_distal=gr_distal.to_dataframe(names=df_distal.columns)
df_tfbs_distal_redundant = df_distal[df_distal["cooperativity"]=="redundant"].reset_index(drop=True)
df_tfbs_distal_codependent = df_distal[df_distal["cooperativity"]=="codependent"].reset_index(drop=True)

df_res_distal_redundant = calc_or_by_threshold(df_tfbs_distal_redundant,trait, snp, n_total_trait, n_total_background)
df_res_distal_codependent = calc_or_by_threshold(df_tfbs_distal_codependent, trait, snp, n_total_trait, n_total_background)


# Create the subplots
fig, axes = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey='row', gridspec_kw={"hspace": 0.4, "wspace": 0.4})




# Plot for proximal redundant (top left)
# error bar no edge
axes[0, 0].errorbar(
    df_res_proximal_redundant["threshold"],
    df_res_proximal_redundant["odds_ratio"],
    yerr=[
        [df_res_proximal_redundant["odds_ratio"][i] - df_res_proximal_redundant["ci_low"][i]
         for i in range(len(df_res_proximal_redundant["odds_ratio"]))],
        [df_res_proximal_redundant["ci_high"][i] - df_res_proximal_redundant["odds_ratio"][i]
         for i in range(len(df_res_proximal_redundant["odds_ratio"]))]
    ],
    fmt='o',
    capsize=0,
    color='tab:blue',
    markersize=5
)
axes[0, 0].set_title("Proximal Redundant")
axes[0, 0].set_ylabel("Odds Ratio")

# Plot for proximal codependent (top right)
axes[0, 1].errorbar(
    df_res_proximal_codependent["threshold"],
    df_res_proximal_codependent["odds_ratio"],
    yerr=[
        [df_res_proximal_codependent["odds_ratio"][i] - df_res_proximal_codependent["ci_low"][i]
         for i in range(len(df_res_proximal_codependent["odds_ratio"]))],
        [df_res_proximal_codependent["ci_high"][i] - df_res_proximal_codependent["odds_ratio"][i]
         for i in range(len(df_res_proximal_codependent["odds_ratio"]))]
    ],
    fmt='o',
    capsize=0,
    color='tab:orange',
    markersize=5
)
axes[0, 1].set_title("Proximal Codependent")

# Plot for distal redundant (bottom left)
axes[1, 0].errorbar(
    df_res_distal_redundant["threshold"],
    df_res_distal_redundant["odds_ratio"],
    yerr=[
        [df_res_distal_redundant["odds_ratio"][i] - df_res_distal_redundant["ci_low"][i]
         for i in range(len(df_res_distal_redundant["odds_ratio"]))],
        [df_res_distal_redundant["ci_high"][i] - df_res_distal_redundant["odds_ratio"][i]
         for i in range(len(df_res_distal_redundant["odds_ratio"]))]
    ],
    fmt='o',
    capsize=0,
    color='tab:blue',
    markersize=5
)
axes[1, 0].set_title("Distal Redundant")
axes[1, 0].set_xlabel("Threshold")
axes[1, 0].set_ylabel("Odds Ratio")

# Plot for distal codependent (bottom right)
axes[1, 1].errorbar(
    df_res_distal_codependent["threshold"],
    df_res_distal_codependent["odds_ratio"],
    yerr=[
        [df_res_distal_codependent["odds_ratio"][i] - df_res_distal_codependent["ci_low"][i]
         for i in range(len(df_res_distal_codependent["odds_ratio"]))],
        [df_res_distal_codependent["ci_high"][i] - df_res_distal_codependent["odds_ratio"][i]
         for i in range(len(df_res_distal_codependent["odds_ratio"]))]
    ],
    fmt='o',
    capsize=0,
    color='tab:orange',
    markersize=5
)
axes[1, 1].set_title("Distal Codependent")
axes[1, 1].set_xlabel("Threshold")

# Adjust layout
plt.tight_layout()
plt.savefig(f"enrichment_by_tf_class_{cell_line}.pdf")
plt.close()
