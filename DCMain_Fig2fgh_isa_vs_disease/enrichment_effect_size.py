import pandas as pd
import numpy as np
from pybedtools import BedTool
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

import sys
sys.path.append("/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import fisher_exact_with_ci

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


prefix="/isdata/alab/people/pcr980/DeepCompare/Pd9_trait_regions"

file_dict={"hepg2":{"clinvar":f"{prefix}/ClinVar_non_benign_non_coding.bed","eqtl":f"{prefix}/liver_eQTLs_PIP_01_non_coding.bed","gwas":f"{prefix}/liver_trait_GWAS_PIP_01_non_coding.bed"},
           "k562":{"clinvar":f"{prefix}/ClinVar_non_benign_non_coding.bed","eqtl":f"{prefix}/whole_blood_eQTLs_PIP_01_non_coding.bed","gwas":f"{prefix}/red_blood_cell_trait_GWAS_PIP_01_non_coding.bed"}
           }
           
track_dict={"hepg2":0,"k562":1}

#--------------------
# Helper functions
#--------------------
def count_overlaps(feature, regions):
    return feature.intersect(regions).to_dataframe().shape[0]



def main(trait, cell_line):

    # Read BED files into BedTool objects
    background_snp = BedTool(f"{prefix}/background_SNPs_non_coding.bed")
    bed_trait = BedTool(file_dict[cell_line][trait])


    # read promoter and enhancer regions
    enhancers = BedTool(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/enhancers_{cell_line}.bed")
    promoters = BedTool(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/promoters_{cell_line}.bed")

    # Read motif info CSV files
    tfbs_enhancers = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_{cell_line}.csv")
    tfbs_promoters = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_{cell_line}.csv")
    df_tfbs= pd.concat([tfbs_enhancers, tfbs_promoters], axis=0, ignore_index=True)
    df_tfbs = df_tfbs.loc[:, ['chromosome', 'start', 'end', 'protein', f'isa_track{track_dict[cell_line]}']]
    # rename column isa_track1 to ISA
    df_tfbs = df_tfbs.rename(columns={"isa_track0": "ISA"})
    df_tfbs = df_tfbs.rename(columns={"isa_track1": "ISA"})
    # select only rows with ISA > 0
    df_tfbs = df_tfbs[df_tfbs["ISA"] > 0].reset_index(drop=True)
    # convert df_tfbs to BedTool


    n_total_background = count_overlaps(enhancers, background_snp) + count_overlaps(promoters, background_snp)
    n_total_trait = count_overlaps(enhancers, bed_trait) + count_overlaps(promoters, bed_trait)


    df_res= []

    for quantile in np.arange(0, 1.0, 0.1):
        threshold = np.quantile(df_tfbs.ISA, quantile)
        df_tfbs_above_thresh = df_tfbs[df_tfbs["ISA"] > threshold].reset_index(drop=True)
        bed_tfbs_above_thresh = BedTool.from_dataframe(df_tfbs_above_thresh)
        n_tfbs_trait = count_overlaps(bed_tfbs_above_thresh, bed_trait)
        n_tfbs_background= count_overlaps(bed_tfbs_above_thresh, background_snp)
        odds_ratio, pval, ci_low, ci_high = fisher_exact_with_ci(np.array([[n_tfbs_trait, n_total_trait - n_tfbs_trait], [n_tfbs_background, n_total_background - n_tfbs_background]]))
        df_res.append({"threshold": quantile, "odds_ratio": odds_ratio, "pvalue": pval, "ci_low": ci_low, "ci_high": ci_high})


    df_res = pd.DataFrame(df_res)

    plt.figure(figsize=(2.1,2.2))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    
    sns.scatterplot(data=df_res, x="threshold", y="odds_ratio", color="orange", s=20)
    plt.errorbar(
        df_res["threshold"],
        df_res["odds_ratio"],
        yerr=[df_res["odds_ratio"] - df_res["ci_low"], df_res["ci_high"] - df_res["odds_ratio"]],
        fmt="none",
        color="orange",
        capsize=0,
        elinewidth=1
    )
    plt.ylim(bottom=0)
    plt.axhline(y=1, color='black', linestyle='--', linewidth=0.5)
    plt.xlabel("Quantile of ISA threshold", fontsize=7)
    plt.ylabel("Odds ratio", fontsize=7)
    plt.xticks(
        ticks=np.arange(0.0, 1.01, 0.1),  # Tick positions from 0.0 to 1.0 in steps of 0.1
        labels=[f"{x:.2f}" for x in np.arange(0.0, 1.01, 0.1)],  # Format as 0.00, 0.10, ..., 1.00
        rotation=90,
        fontsize=5
    )
    plt.yticks(fontsize=5)

    # Ensure tight layout and save the figure
    plt.title(f"{trait} {cell_line}", fontsize=7)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(f"{trait}_{cell_line}.pdf")
    plt.close()


for trait in ["clinvar", "eqtl", "gwas"]:
    for cell_line in ["hepg2", "k562"]:
        main(trait, cell_line)
        logger.info(f"Finished {trait} {cell_line}")
        
        
# nohup python3 enrichment_effect_size.py > enrichment_effect_size.log &