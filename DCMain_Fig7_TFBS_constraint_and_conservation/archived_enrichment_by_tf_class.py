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

dict_track_num={"hepg2":0, "k562":1}


cds=BedTool("cds.bed")
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




def calc_or(df_tfbs, bed_trait, snp, n_total_trait, n_total_background):
    df_tfbs_above_thresh = df_tfbs[df_tfbs[f"isa_track{dict_track_num[cell_line]}"] > 0].reset_index(drop=True)
    bed_tfbs_above_thresh = BedTool.from_dataframe(df_tfbs_above_thresh)
    n_tfbs_trait = count_overlaps(bed_tfbs_above_thresh, bed_trait)
    n_tfbs_background= count_overlaps(bed_tfbs_above_thresh, snp)
    odds_ratio, pval, ci_low, ci_high = fisher_exact_with_ci(np.array([[n_tfbs_trait, n_total_trait - n_tfbs_trait], [n_tfbs_background, n_total_background - n_tfbs_background]]))
    return odds_ratio, pval, ci_low, ci_high




def get_all_dhs(cell_line):
    dhs_proximal=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/dhs_proximal_{cell_line}.tsv")
    dhs_distal=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/dhs_distal_{cell_line}.tsv")
    dhs=pd.concat([dhs_proximal, dhs_distal], axis=0, ignore_index=True)
    # convert to bedtool
    dhs=BedTool.from_dataframe(dhs)
    return dhs






def preprocess(file_name):
    prefix_motif_info="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs"
    df_tfbs = pd.read_csv(f"{prefix_motif_info}_{file_name}.csv")
    df_tfbs = df_tfbs.drop(columns=["Unnamed: 0"])
    if "hepg2" in file_name:
        cell_line="hepg2"
    elif "k562" in file_name:
        cell_line="k562"
    # select only constrained TFs (z>0)
    df_constraint=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_{cell_line}_dhs.csv")
    tfs_constrained=df_constraint[df_constraint["z"]>0]["protein"].tolist()
    df_tfbs=df_tfbs[df_tfbs["protein"].isin(tfs_constrained)].reset_index(drop=True)
    # add cooperativity info
    df_tfbs=add_tf_codependency(df_tfbs, cell_line)
    return df_tfbs







for cell_line in ["hepg2","k562"]:        
    # load tfbs data related to cell line
    df_tfbs_proximal = preprocess(f"proximal_{cell_line}")
    df_tfbs_distal = preprocess(f"distal_{cell_line}")
    #
    df_tfbs_proximal_redundant = df_tfbs_proximal[df_tfbs_proximal["cooperativity"]=="redundant"].reset_index(drop=True)
    df_tfbs_proximal_codependent = df_tfbs_proximal[df_tfbs_proximal["cooperativity"]=="codependent"].reset_index(drop=True)
    #
    df_tfbs_distal_redundant = df_tfbs_distal[df_tfbs_distal["cooperativity"]=="redundant"].reset_index(drop=True)
    df_tfbs_distal_codependent = df_tfbs_distal[df_tfbs_distal["cooperativity"]=="codependent"].reset_index(drop=True)
    # Calculate overlaps as base
    dhs = get_all_dhs(cell_line)
    n_total_background = count_overlaps(snp, dhs)
    df_res=pd.DataFrame()
    for trait_indicator in ["eqtl","clinvar","gwas"]:
        trait = BedTool(file_dict[cell_line][trait_indicator])
        n_total_trait = count_overlaps(trait, dhs)
        #
        res_proximal_redundant= calc_or(df_tfbs_proximal_redundant, trait, snp, n_total_trait, n_total_background)
        res_proximal_codependent = calc_or(df_tfbs_proximal_codependent, trait, snp, n_total_trait, n_total_background)
        #
        res_distal_redundant = calc_or(df_tfbs_distal_redundant,trait, snp, n_total_trait, n_total_background)
        res_distal_codependent = calc_or(df_tfbs_distal_codependent, trait, snp, n_total_trait, n_total_background)
        # write results to a df
        df_new=pd.DataFrame({"proximal_redundant":res_proximal_redundant, "proximal_codependent":res_proximal_codependent, "distal_redundant":res_distal_redundant, "distal_codependent":res_distal_codependent}, index=["odds_ratio", "pval", "ci_low", "ci_high"])
        df_new["trait"]=trait_indicator
        df_res=pd.concat([df_res, df_new], axis=0)
    df_res.to_csv(f"enrichment_by_tf_class_{cell_line}.csv")





# nohup python3 enrichment_by_tf_class.py &