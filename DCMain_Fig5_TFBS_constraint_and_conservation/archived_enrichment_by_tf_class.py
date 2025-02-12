import pandas as pd
import numpy as np
from pybedtools import BedTool
from matplotlib import pyplot as plt
import seaborn as sns
import warnings


import sys
sys.path.append("/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import fisher_exact_with_ci


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42

prefix="/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions"

file_dict={"hepg2":{"clinvar":f"{prefix}/ClinVar_non_benign_non_coding.bed","eqtl":f"{prefix}/liver_eQTLs_PIP_01_non_coding.bed","gwas":f"{prefix}/liver_trait_GWAS_PIP_01_non_coding.bed"},
           "k562":{"clinvar":f"{prefix}/ClinVar_non_benign_non_coding.bed","eqtl":f"{prefix}/whole_blood_eQTLs_PIP_01_non_coding.bed","gwas":f"{prefix}/red_blood_cell_trait_GWAS_PIP_01_non_coding.bed"}
           }

dict_track_num={"hepg2":6, "k562":7}


cds=BedTool("cds.bed")
# change to df
snp = BedTool("/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions/background_SNPs_non_coding.bed")





def add_tf_codependency(df,suffix):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{suffix}_dhs.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{suffix}_dhs.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df



def count_overlaps(feature, regions):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")  
        return feature.intersect(regions).to_dataframe().shape[0]




def calc_or_by_threshold(df_tfbs, bed_trait, snp, n_total_trait, n_total_background):
    df_res = []
    for threshold in np.arange(0, 0.4, 0.05):
        df_tfbs_above_thresh = df_tfbs[df_tfbs[f"isa_track{dict_track_num[cell_line]}"] > threshold].reset_index(drop=True)
        bed_tfbs_above_thresh = BedTool.from_dataframe(df_tfbs_above_thresh)
        n_tfbs_trait = count_overlaps(bed_tfbs_above_thresh, bed_trait)
        n_tfbs_background= count_overlaps(bed_tfbs_above_thresh, snp)
        print(n_tfbs_trait, n_tfbs_background)
        odds_ratio, pval, ci_low, ci_high = fisher_exact_with_ci(np.array([[n_tfbs_trait, n_total_trait - n_tfbs_trait], [n_tfbs_background, n_total_background - n_tfbs_background]]))
        df_res.append({"threshold": threshold, "odds_ratio": odds_ratio, "pvalue": pval, "ci_low": ci_low, "ci_high": ci_high})
    return pd.DataFrame(df_res)




def get_all_dhs(cell_line):
    dhs_proximal=pd.read_csv(f"Pd1_dnase_annotated/dhs_proximal_{cell_line}.csv")
    dhs_distal=pd.read_csv(f"Pd1_dnase_annotated/dhs_distal_{cell_line}.csv")
    dhs=pd.concat([dhs_proximal, dhs_distal], axis=0, ignore_index=True)
    # convert to bedtool
    dhs=BedTool.from_dataframe(dhs)
    return dhs





def preprocess(df,cell_line):
    df=df.iloc[:,1:]
    df= add_tf_codependency(df,f"{cell_line}")
    gr=BedTool.from_dataframe(df)
    gr=gr.intersect(cds, v=True)
    df=gr.to_dataframe(names=df.columns)
    return df




def plot(df_proximal_redundant,df_proximal_codependent,df_distal_redundant,df_distal_codependent,cell_line,trait_indicator):
    df_proximal_redundant["category"] = "Proximal Redundant"
    df_proximal_codependent["category"] = "Proximal Codependent"
    df_distal_redundant["category"] = "Distal Redundant"
    df_distal_codependent["category"] = "Distal Codependent"
    # Combine all dataframes into one
    plot_data = pd.concat([
        df_proximal_redundant,
        df_proximal_codependent,
        df_distal_redundant,
        df_distal_codependent
    ], ignore_index=True)
    category_info = {
        "Proximal Redundant": (0, 0, "tab:blue"),
        "Proximal Codependent": (0, 1, "tab:orange"),
        "Distal Redundant": (1, 0, "tab:blue"),
        "Distal Codependent": (1, 1, "tab:orange"),
    }
    fig, axes = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey='row')
    # Plot for proximal redundant (top left)
    for category, group_df in plot_data.groupby("category"):
        row, col, color = category_info[category]
        axes[row, col].errorbar(
            group_df["threshold"],
            group_df["odds_ratio"],
            yerr=[
                group_df["odds_ratio"] - group_df["ci_low"],
                group_df["ci_high"] - group_df["odds_ratio"]
            ],
            fmt='o',
            capsize=0,
            color=color,
            markersize=5
        )
        axes[row, col].set_title(category)
        if row == 1:
            axes[row, col].set_xlabel("Threshold")
        if col == 0:
            axes[row, col].set_ylabel("Odds Ratio")
    # Adjust layout
    plt.suptitle(f"Enrichment by TF Class ({cell_line} {trait_indicator})")
    plt.tight_layout()
    plt.savefig(f"enrichment_by_tf_class_sure_{cell_line}_{trait_indicator}.pdf")
    plt.close()



prefix_motif_info="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs"
for cell_line in ["hepg2","k562"]:
    # load tfbs data related to cell line
    df_tfbs_proximal = pd.read_csv(f"{prefix_motif_info}_proximal_{cell_line}.csv")
    df_tfbs_distal = pd.read_csv(f"{prefix_motif_info}_distal_{cell_line}.csv")
    # preprocess
    df_tfbs_proximal = preprocess(df_tfbs_proximal,cell_line)
    df_tfbs_distal = preprocess(df_tfbs_distal,cell_line)
    #
    df_tfbs_proximal_redundant = df_tfbs_proximal[df_tfbs_proximal["cooperativity"]=="redundant"].reset_index(drop=True)
    df_tfbs_proximal_codependent = df_tfbs_proximal[df_tfbs_proximal["cooperativity"]=="codependent"].reset_index(drop=True)
    #
    df_tfbs_distal_redundant = df_tfbs_distal[df_tfbs_distal["cooperativity"]=="redundant"].reset_index(drop=True)
    df_tfbs_distal_codependent = df_tfbs_distal[df_tfbs_distal["cooperativity"]=="codependent"].reset_index(drop=True)
    # Calculate overlaps as base
    dhs = get_all_dhs(cell_line)
    n_total_background = count_overlaps(snp, dhs)
    #
    for trait_indicator in ["eqtl","clinvar","gwas"]:
        trait = BedTool(file_dict[cell_line][trait_indicator])
        n_total_trait = count_overlaps(trait, dhs)
        #
        df_res_proximal_redundant = calc_or_by_threshold(df_tfbs_proximal_redundant, trait, snp, n_total_trait, n_total_background)
        df_res_proximal_codependent = calc_or_by_threshold(df_tfbs_proximal_codependent, trait, snp, n_total_trait, n_total_background)
        #
        df_res_distal_redundant = calc_or_by_threshold(df_tfbs_distal_redundant,trait, snp, n_total_trait, n_total_background)
        df_res_distal_codependent = calc_or_by_threshold(df_tfbs_distal_codependent, trait, snp, n_total_trait, n_total_background)
        #
        #
        plot(df_res_proximal_redundant,df_res_proximal_codependent,df_res_distal_redundant,df_res_distal_codependent,cell_line,trait_indicator)




# nohup python3 enrichment_by_tf_class.py &