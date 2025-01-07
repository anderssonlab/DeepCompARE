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
        warnings.filterwarnings("ignore")  # Suppress warnings within this block
        print(feature.intersect(regions).to_dataframe())
        return feature.intersect(regions).to_dataframe().shape[0]



cell_line="k562"
dict_track_num={"hepg2":0, "k562":1}

# Read BED files into BedTool objects
cds=BedTool("cds.bed")
eqtl = BedTool("/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions/whole_blood_eQTLs_PIP_01_non_coding.bed")
snp = BedTool("/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions/background_SNPs_non_coding.bed")


dhs_proximal=pd.read_csv(f"Pd1_dnase_annotated/dhs_proximal_{cell_line}.csv")
dhs_distal=pd.read_csv(f"Pd1_dnase_annotated/dhs_distal_{cell_line}.csv")
dhs=pd.concat([dhs_proximal, dhs_distal], axis=0, ignore_index=True)
dhs = BedTool.from_dataframe(dhs)


# Load TFBS DHS data
df_proximal = pd.read_csv(f"motif_info_thresh_500_dhs_proximal_{cell_line}.csv")
df_distal = pd.read_csv(f"motif_info_thresh_500_dhs_distal_{cell_line}.csv")
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
n_total_trait = count_overlaps(eqtl, dhs)



df_tfbs = gr_proximal.to_dataframe(names=df_proximal.columns)
df_tfbs=df_tfbs[df_tfbs["cooperativity"]=="redundant"].reset_index(drop=True)


def calc_or_by_threshold(df_tfbs, bed_trait, snp, n_total_trait, n_total_background):
    df_res = []
    for threshold in np.arange(0, 0.5, 0.02):
        df_tfbs_above_thresh = df_tfbs[df_tfbs[f"isa_track{dict_track_num[cell_line]}"] > threshold].reset_index(drop=True)
        print(df_tfbs_above_thresh.shape)
        bed_tfbs_above_thresh = BedTool.from_dataframe(df_tfbs_above_thresh)
        n_tfbs_trait = count_overlaps(bed_tfbs_above_thresh, bed_trait)
        print(n_tfbs_trait)
        n_tfbs_background= count_overlaps(bed_tfbs_above_thresh, snp)
        odds_ratio, pval, ci_low, ci_high = fisher_exact_with_ci(np.array([[n_tfbs_trait, n_total_trait - n_tfbs_trait], [n_tfbs_background, n_total_background - n_tfbs_background]]))
        df_res.append({"threshold": threshold, "odds_ratio": odds_ratio, "pvalue": pval, "ci_low": ci_low, "ci_high": ci_high})
    return pd.DataFrame(df_res)




df_proximal = gr_proximal.to_dataframe(names=df_proximal.columns)
df_tfbs_proximal_redundant = df_proximal[df_proximal["cooperativity"]=="redundant"].reset_index(drop=True)
df_tfbs_proximal_codependent = df_proximal[df_proximal["cooperativity"]=="codependent"].reset_index(drop=True)
df_tfbs_proximal_unknown = df_proximal[df_proximal["cooperativity"]=="unknown"].reset_index(drop=True)

df_res_proximal_redundant = calc_or_by_threshold(df_tfbs_proximal_redundant, eqtl, snp, n_total_trait, n_total_background)
df_res_proximal_codependent = calc_or_by_threshold(df_tfbs_proximal_codependent, eqtl, snp, n_total_trait, n_total_background)
df_res_proximal_unknown = calc_or_by_threshold(df_tfbs_proximal_unknown, eqtl, snp, n_total_trait, n_total_background)





df_distal=gr_distal.to_dataframe(names=df_distal.columns)
df_tfbs_distal_redundant = df_distal[df_distal["cooperativity"]=="redundant"].reset_index(drop=True)
df_tfbs_distal_codependent = df_distal[df_distal["cooperativity"]=="codependent"].reset_index(drop=True)
df_tfbs_distal_unknown = df_distal[df_distal["cooperativity"]=="unknown"].reset_index(drop=True)


df_res_distal_redundant = calc_or_by_threshold(df_tfbs_distal_redundant,eqtl, snp, n_total_trait, n_total_background)
df_res_distal_codependent = calc_or_by_threshold(df_tfbs_distal_codependent, eqtl, snp, n_total_trait, n_total_background)
df_res_distal_unknown = calc_or_by_threshold(df_tfbs_distal_unknown, eqtl, snp, n_total_trait, n_total_background)



# Visualization
plt.figure(figsize=(10, 6))

sns.barplot(x="Locus", y="eQTL_P", hue="Type", data=results_df)
plt.axhline(1, linestyle="--", color="gray", linewidth=0.5)
plt.xlabel("ISA Class")
plt.ylabel("Odds Ratio")
plt.title("Enrichment vs ISA DHS Cage Location Source")
plt.show()
