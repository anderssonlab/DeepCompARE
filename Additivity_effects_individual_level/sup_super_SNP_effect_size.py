import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from loguru import logger
from scipy import stats


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_and_label, calc_or, transform_data_for_plotting, plot_or


#-----------------------------------------------------
# helper functions
#-----------------------------------------------------
def get_contextual_tf_type(df,dataset):
    df_tf_property=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_additivity_property.csv")
    sub_tfs = df_tf_property[df_tf_property[dataset]=="sub"]["protein"].to_list()
    super_tfs = df_tf_property[df_tf_property[dataset]=="super"]["protein"].to_list()
    tf_type= np.where(df["motif"].isin(sub_tfs),"sub",np.where(df["motif"].isin(super_tfs),"super","other"))
    tf_type=pd.Categorical(tf_type,categories=["sub","super","other"])
    return tf_type

def read_tfbs_maf(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_{file_suffix}.csv",header=None)
    df.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]
    df["dataset"]=file_suffix
    df["log10_AF"]=np.log10(df["AF"])
    return df

def read_maf_effect_size(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Effect_size_vs_maf/maf_with_effect_size_{file_suffix}.csv",header=None,index_col=0)
    df.reset_index(drop=True,inplace=True)
    df.columns=["Chromosome","Start","End","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
    return df

def join_tfbs_maf_effect_size(file_suffix):
    tfbs_maf=read_tfbs_maf(file_suffix)
    maf_effect_size=read_maf_effect_size(file_suffix)
    df=tfbs_maf.join(maf_effect_size.set_index(["Chromosome","Start","End","ID","REF","ALT","AF"]),on=["Chromosome","Start","End","ID","REF","ALT","AF"],how="inner")
    return df

#-----------------------------------------------------
# read in all SNPs, allele frequency, TFBS and effect size
#-----------------------------------------------------


# df_tfbs_maf_es is for effect size analysis
df_tfbs_maf_es=pd.concat([join_tfbs_maf_effect_size("enhancers_hepg2"),
                          join_tfbs_maf_effect_size("promoters_hepg2"),
                          join_tfbs_maf_effect_size("enhancers_k562"),
                          join_tfbs_maf_effect_size("promoters_k562")],
                          axis=0).reset_index(drop=True)
df_tfbs_maf_es=df_tfbs_maf_es[df_tfbs_maf_es["chip_evidence"]==True].reset_index(drop=True)
df_tfbs_maf_es=bin_and_label(df_tfbs_maf_es, "log10_AF", [-np.inf,-4, -2, 0], "AF_bin")




#--------------------------------------------------------------------------------------------
# Allele frequency analysis, sub v.s. rest
#--------------------------------------------------------------------------------------------
sub_tfs = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/sub_tfs.txt",header=None)[0].to_list()
df_tfbs_maf_es["TF_type"]=np.where(df_tfbs_maf_es["motif"].isin(sub_tfs),"sub","other")
# df_tfbs_maf_es["TF_type"]=np.where(df_tfbs_maf_es["motif"].isin(sub_tfs),"sub",np.where(df_tfbs_maf_es["motif"].isin(super_tfs),"super","other"))

# plot distribution of log10_AF, hue by TF_type
sns.kdeplot(data=df_tfbs_maf_es,x="log10_AF",hue="TF_type",common_norm=False)
plt.title("Allele frequency distribution of SNVs overlapping REs")
plt.savefig("Plots_maf/allele_frequency_distribution.pdf")

sns.histplot(data=df_tfbs_maf_es,x="log10_AF",hue="TF_type",bins=100)
plt.title(f"Allele frequency of SNVs overlapping REs of all")
plt.yscale("log")
plt.savefig(f"Plots_maf/allele_frequency_distribution_all.pdf")
plt.close()
# plot percentage
sns.histplot(data=df_tfbs_maf_es,x="log10_AF",hue="TF_type",bins=100,multiple="fill")
plt.title(f"Allele frequency of SNVs overlapping REs of all")
plt.savefig(f"Plots_maf/allele_frequency_distribution_fill_all.pdf")
plt.close()








# for 4 datasets combined
df_es_binned = df_tfbs_maf_es.loc[:, ["AF_bin", "TF_type"]].copy().groupby(["AF_bin", "TF_type"]).size().unstack()
df_es_binned = calc_or(df_es_binned)
plot_or(transform_data_for_plotting(df_es_binned, "AF_bin"),
        "AF_bin", "odds_ratio", "4 Files combined",{'sub': "#1f77b4", 'super': '#ff7f0e'},
        out_name="Plots_maf/af_enrichment_all.pdf",
        ax=None)





# for each dataset
for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
    df_subset=df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
    df_subset["TF_type"]=get_tf_type(df_subset,dataset)
    df_es_binned = df_subset.loc[:, ["AF_bin", "TF_type"]].copy().groupby(["AF_bin", "TF_type"]).size().unstack()
    df_es_binned = calc_or(df_es_binned)
    plot_or(transform_data_for_plotting(df_es_binned, "AF_bin"),
            "AF_bin", "odds_ratio",dataset,{'sub': "#1f77b4", 'super': '#ff7f0e'},
            out_name=f"Plots_maf/af_enrichment_{dataset}.pdf",
            ax=None)


for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
    df_subset=df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
    df_subset["TF_type"]=get_tf_type(df_subset,dataset)
    # plot histogram
    sns.histplot(data=df_subset,x="log10_AF",hue="TF_type",bins=100)
    plt.title(f"Allele frequency of SNVs overlapping REs of {dataset}")
    plt.yscale("log")
    plt.savefig(f"Plots_maf/allele_frequency_distribution_{dataset}.pdf")
    plt.close()
    # plot percentage
    sns.histplot(data=df_subset,x="log10_AF",hue="TF_type",bins=100,multiple="fill")
    plt.title(f"Allele frequency of SNVs overlapping REs of {dataset}")
    plt.savefig(f"Plots_maf/allele_frequency_distribution_fill_{dataset}.pdf")
    plt.close()



#--------------------------------------------------------------------------------------------
# enrichment of sub and super
#--------------------------------------------------------------------------------------------




def plot_combined_es_enrichment(dataset, df_plot, x_colname, y_colname, track_num, color_mapping):
    af_bins = df_plot['AF_bin'].unique()
    fig, axs = plt.subplots(len(af_bins), 1, figsize=(8, 8), sharex=True)
    for i, af_bin in enumerate(af_bins):
        df_es_binned = df_plot[df_plot["AF_bin"] == af_bin].reset_index(drop=True)
        df_es_binned = df_es_binned.loc[:, ["Bin", "TF_type"]].copy().groupby(["Bin", "TF_type"]).size().unstack()
        df_es_binned = calc_or(df_es_binned)
        df = transform_data_for_plotting(df_es_binned, "Bin")
        plot_or(df, x_colname, y_colname, f"log10 AF: {af_bin}", color_mapping, ax=axs[i])
    axs[-1].set_xlabel(x_colname)
    fig.text(0.04, 0.5, y_colname, va='center', rotation='vertical')
    plt.suptitle(f"{dataset}, track {track_num}")
    plt.tight_layout()
    plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.15)
    plt.savefig(f"Plots_maf/es_enrichment_{dataset}_track{track_num}.pdf")
    plt.close()
    


color_mapping = {'sub': "#1f77b4", 'super': '#ff7f0e'}


for dataset in ["enhancers_hepg2","promoters_hepg2","enhancers_k562","promoters_k562"]:
    df_subset = df_tfbs_maf_es[df_tfbs_maf_es["dataset"]==dataset].reset_index(drop=True)
    df_subset["TF_type"]=get_tf_type(df_subset,dataset)
    if "hepg2" in dataset:
        track_list = [0, 2, 4, 6]
    else:
        track_list = [1, 3, 5, 7]
    for track_num in track_list:
        df_plot = bin_and_label(df_subset, f"track_{track_num}", [-np.inf, -0.2, -0.1, 0, 0.1, 0.2, np.inf])
        plot_combined_es_enrichment(dataset, df_plot, "Bin", "odds_ratio", track_num, color_mapping)
















for cell_type in ["hepg2","k562"]:
    df_subset=df_tfbs_maf_es[df_tfbs_maf_es["dataset"].str.contains(cell_type)].reset_index(drop=True)
    if cell_type=="hepg2":
        track_list=[0,2,4,6]
    else:
        track_list=[1,3,5,7]
    for track_num in track_list:
        df_subset=bin_and_label(df_subset,f"track_{track_num}", [-np.inf, -0.2, -0.1, 0, 0.1, 0.2, np.inf])
        for i, af_bin in enumerate(df_es_binned['AF_bin'].unique()):
            logger.info(f"Processing track {track_num}, allele frequency strata: {af_bin}")
            df_es_binned=df_subset[df_subset["AF_bin"]==af_bin].reset_index(drop=True)
            df_es_binned=df_es_binned.loc[:,["Bin","TF_type"]].copy().groupby(["Bin","TF_type"]).size().unstack()
            df_es_binned=calc_or(df_es_binned)
            df=transform_data_for_plotting(df_es_binned,"Bin") 
            plot_or(df,"Bin","odds_ratio",
                    f"Allele frequency strata: {af_bin}, track {track_num}",
                    f"Plots_maf/es_enrichment_{af_bin}_track{track_num}.pdf",
                    {'sub': "#1f77b4", 'super': '#ff7f0e'})




















#---------------------------------------
# Archived
#---------------------------------------



# def get_additivity_profile(df):
#     df_additivity=pd.read_csv("tf_property_long_format.csv")
#     df_additivity.rename(columns={"protein":"motif"},inplace=True)
#     df=df.merge(df_additivity,on=["dataset","motif"],how="inner")
#     df=df[df["tf_property"].isin(["super","sub"])].reset_index(drop=True)
#     df["tf_property"]=pd.Categorical(df["tf_property"],categories=["sub","super"])
#     return df

# def plot_contextual_sub_vs_super(df_orig, value_col):
#     df=df_orig.copy()
#     df[value_col]=df[value_col].abs()
#     for dataset in df["dataset"].unique():
#         df_subset = df[df["dataset"] == dataset].reset_index(drop=True).copy()
#         _, p_value = stats.mannwhitneyu(df_subset[df_subset["tf_property"] == "sub"][value_col],
#                                              df_subset[df_subset["tf_property"] == "super"][value_col])
#         sub_median = df_subset[df_subset["tf_property"] == "sub"][value_col].median()
#         super_median = df_subset[df_subset["tf_property"] == "super"][value_col].median()
#         larger_median_group = "sub" if sub_median > super_median else "super"
#         sns.kdeplot(data=df_subset, x=value_col, hue="tf_property", common_norm=False)
#         plt.title(f"{value_col} {dataset}")
#         plt.annotate(f"p={p_value:.1e}", xy=(0.5, 0.5), xycoords="axes fraction")
#         plt.annotate(f"larger median: {larger_median_group}", xy=(0.5, 0.4), xycoords="axes fraction")
#         plt.savefig(f"Plots_maf/contextual_sub_vs_super_{value_col}_{dataset}.pdf")
#         plt.close()


# df_tfbs_maf_es=get_additivity_profile(df_tfbs_maf_es)



# plot_contextual_sub_vs_super(df_tfbs_maf, "log10_AF")
# for i in range(8):
#     plot_contextual_sub_vs_super(df_tfbs_maf_es, f"track_{i}")


