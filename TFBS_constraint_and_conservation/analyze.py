import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu,pearsonr
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_two_columns, calc_or, plot_or
from utils import get_track_num
from seq_ops import SeqExtractor



#-------------------
# Helper functions
#-------------------

def read_file(file_suffix,seq_extractor):
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv')
    track_nums=get_track_num(file_suffix)
    track_nums_remove=[i for i in range(16) if i not in track_nums]
    cols_remove=[f"ism_track{i}" for i in track_nums_remove]
    df.drop(cols_remove, axis=1, inplace=True)
    mapper={0:"cage",1:"cage",2:"dhs",3:"dhs",4:"starr",5:"starr",6:"sure",7:"sure"}
    df.rename(columns={f"ism_track{i}":f"ism_{mapper[i]}" for i in track_nums}, inplace=True)
    cols_remove=[col for col in df.columns if col.startswith('pred_orig')]
    df.drop(cols_remove, axis=1, inplace=True)
    df["motif"]=df.apply(lambda row: seq_extractor.get_seq((row["chromosome"],row["start"],row["end"])), axis=1)
    df["gc_content"]=df["motif"].apply(lambda x: (x.count("G")+x.count("C"))/len(x))
    df["max_af"] = df["af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["min_af"] = df["af"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["num_variants"] = df["af"].apply(lambda x: len(str(x).split(":")))
    df["num_common_variants"] = df["af"].apply(lambda x: sum([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["241way_max"] = df["241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    return df



def add_tf_cooperativity(df):
    tfs_codependent=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_codependent.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_redundant.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df



seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")



#--------------------------------------------------------------
# Analysis 1: effect size v.s. allele frequency
# Goal: reproduce SNP level enrichment
#--------------------------------------------------------------


# Read in data
df_promoter_hepg2=read_file("promoters_hepg2",seq_extractor)
df_promoter_hepg2=add_tf_cooperativity(df_promoter_hepg2)
df_promoter_hepg2["dataset"]="promoters_hepg2"
df_enhancer_hepg2=read_file("enhancers_hepg2",seq_extractor)
df_enhancer_hepg2=add_tf_cooperativity(df_enhancer_hepg2)
df_enhancer_hepg2["dataset"]="enhancers_hepg2"
df_promoter_k562=read_file("promoters_k562",seq_extractor)
df_promoter_k562=add_tf_cooperativity(df_promoter_k562)
df_promoter_k562["dataset"]="promoters_k562"
df_enhancer_k562=read_file("enhancers_k562",seq_extractor)
df_enhancer_k562=add_tf_cooperativity(df_enhancer_k562)
df_enhancer_k562["dataset"]="enhancers_k562"


df=pd.concat([df_promoter_hepg2,df_enhancer_hepg2,df_promoter_k562,df_enhancer_k562], axis=0)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])

# TODO: maybe make more bins
# plot effect size v.s. max allele frequency
def plot_or_max_af(df,dataset,track_name,suffix=None):
    df_plot=bin_two_columns(df,f"ism_{track_name}",[0, 0.1, 0.2, 0.5, np.inf],
                       "max_af",{"0 - 0.001": "max(AF)<0.001", "0.001 - 0.01": "low", "0.01 - 1":"max(AF)>0.01"},
                       [0,0.001,0.01,1],"TFBS_type")
    df_plot=df_plot[df_plot.sum(axis=1)>10]
    df_plot=calc_or(df_plot,"ISA score","TFBS_type",out_group="low")
    plot_or(df_plot, 'ISA score', 'odds_ratio', "TFBS_type",
            f"Odds ratio ({dataset}, track {track_name})",
            {'max(AF)<0.001': "#1f77b4", 'max(AF)>0.01': '#ff7f0e'},
            f"Plots/or_max_af_{dataset}_{track_name}.pdf")


for dataset in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
    for track_name in ["cage","dhs","starr","sure"]:
        logger.info(f"Processing {dataset}, track {track_name}")
        df_sub=df[df["dataset"]==dataset].reset_index(drop=True)
        df_sub=df_sub[df_sub[f"ism_{track_name}"]>0].reset_index(drop=True)
        plot_or_max_af(df_sub,dataset,track_name)




# plot effect size v.s. max phylop score
def plot_or_241way_max(df,dataset,track_name,suffix=None):
    df_plot=bin_two_columns(df,f"ism_{track_name}",[0, 0.1, 0.2, 0.5, np.inf],
                       "241way_max",{"0 - 1": "max(phyloP)<1", "1 - 3": "conserved", "3 - inf":"max(phyloP)>3"},
                       [0,1,3,np.inf],"TFBS_type")
    df_plot=df_plot[df_plot.sum(axis=1)>10]
    df_plot=calc_or(df_plot,"ISA score","TFBS_type",out_group="low")
    plot_or(df_plot, 'ISA score', 'odds_ratio', "TFBS_type",
            f"Odds ratio ({dataset}, track {track_name})",
            {'max(phyloP)<1': "#1f77b4", 'max(phyloP)>3': '#ff7f0e'},
            f"Plots/or_241way_max_{dataset}_{track_name}.pdf")



for dataset in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
    for track_name in ["cage","dhs","starr","sure"]:
        logger.info(f"Processing {dataset}, track {track_name}")
        df_sub=df[df["dataset"]==dataset].reset_index(drop=True)
        df_sub=df_sub[df_sub[f"ism_{track_name}"]>0].reset_index(drop=True)
        plot_or_241way_max(df_sub,dataset,track_name)




#--------------------------------------------------------------
# Analysis 2: merge promoters and enhancers, compare cooperativity
#--------------------------------------------------------------

def preprocess_by_cell_type(df):
    df=df.groupby("protein").agg({"241way_max":"mean",
                                  "241way_min":"mean",
                                  "gc_content":"mean",
                                  "num_expected_variants":"sum",
                                  "num_rare_variants":"sum",
                                }).reset_index()
    df=add_tf_cooperativity(df)
    return df


cell="k562"
df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
df_sub=preprocess_by_cell_type(df_sub)
df_sub=df_sub[df_sub["cooperativity"]!="unknown"].reset_index(drop=True)

g = sns.jointplot(data=df_sub, y="241way_max", x="z", hue="cooperativity", marginal_kws={'common_norm': False})
plt.subplots_adjust(top=0.9)  # Adjust the top to make space for the title
g.fig.suptitle(f'{cell}', fontsize=16)
plt.savefig(f"Plots/joint_241way_max_z_{cell}.pdf")
plt.close()

#--------------------------------------------------------------
# Analysis 3: separate by cooperativity, compare promoters and enhancers
#--------------------------------------------------------------


def preprocess_group_by_tf(file_suffix,seq_extractor):
    df=read_file(file_suffix,seq_extractor)
    df=df.groupby("protein").agg({"have_nonrare_variant":"mean",
                                  "241way_max":"mean",
                                  "241way_min":"mean",
                                  "gc_content":"mean",
                                  "num_expected_variants":"sum",
                                  "num_rare_variants":"sum",
                                  "num_variants":"sum",
                                  "num_common_variants":"sum"
                                }).reset_index()
    df=add_tf_cooperativity(df)
    df=calc_gnocchi(df)
    return df



df_promoter_hepg2=preprocess_group_by_tf("promoters_hepg2",seq_extractor)
df_promoter_hepg2["dataset"]="promoters_hepg2"
df_enhancer_hepg2=preprocess_group_by_tf("enhancers_hepg2",seq_extractor)
df_enhancer_hepg2["dataset"]="enhancers_hepg2"
df_promoter_k562=preprocess_group_by_tf("promoters_k562",seq_extractor)
df_promoter_k562["dataset"]="promoters_k562"
df_enhancer_k562=preprocess_group_by_tf("enhancers_k562",seq_extractor)
df_enhancer_k562["dataset"]="enhancers_k562"


df=pd.concat([df_promoter_hepg2,df_enhancer_hepg2,df_promoter_k562,df_enhancer_k562], axis=0)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])


xcol="z"
ycol="241way_max"
for cell in ["hepg2","k562"]:
    for coop in ["codependent","redundant","all"]:
        logger.info(f"Processing {cell} {coop}")
        df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
        if coop!="all":
            df_sub=df_sub[df_sub["cooperativity"]==coop].reset_index(drop=True)
        ustat, pval = mannwhitneyu(df_sub[df_sub["re"]=="promoters"][xcol], df_sub[df_sub["re"]=="enhancers"][xcol])
        logger.info(f"U-statistic: {ustat}, p-value: {pval}")
        logger.info(f"Mean z score: {df_sub[df_sub['re']=='promoters'][xcol].mean()} {df_sub[df_sub['re']=='enhancers'][xcol].mean()}")
        g = sns.jointplot(data=df_sub, x=xcol, y=ycol,hue="re", marginal_kws={'common_norm': False})
        plt.subplots_adjust(top=0.9)  # Adjust the top to make space for the title
        g.fig.suptitle(f'{coop} tfs({cell})', fontsize=16)
        plt.savefig(f"Plots/joint_{ycol}_{xcol}_promoters_vs_enhancers_{coop}_{cell}.pdf")
        plt.close()






col="z"
cell="k562"
logger.info(f"Processing {cell} {col}")
df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
df_redundant=df_sub[df_sub["cooperativity"]=="redundant"].reset_index(drop=True)
df_codependent=df_sub[df_sub["cooperativity"]=="codependent"].reset_index(drop=True)
# pivot to compare z between promoters and enhancers
df_redundant=df_redundant.pivot(index="protein", columns="re", values=col).dropna().reset_index()
df_codependent=df_codependent.pivot(index="protein", columns="re", values=col).dropna().reset_index()
# how many enhancers>promoters
logger.info(f"Number of proteins where enhancers>promoters (redundant tfs): {np.mean(df_redundant['enhancers']>df_redundant['promoters'])}")
logger.info(f"Number of proteins where enhancers>promoters (codependent tfs): {np.mean(df_codependent['enhancers']>df_codependent['promoters'])}")
df_redundant["cooperativity"]="redundant"
df_codependent["cooperativity"]="codependent"
df_pivot=pd.concat([df_redundant,df_codependent], axis=0)

plt.figure(figsize=(6,6))
sns.scatterplot(x="promoters", y="enhancers", data=df_pivot, hue="cooperativity")
plt.title(f'{col}, {cell}')
# abline
min_val=min(df_pivot["promoters"].min(),df_pivot["enhancers"].min())
max_val=max(df_pivot["promoters"].max(),df_pivot["enhancers"].max())
plt.plot([min_val,max_val],[min_val,max_val], color="black")
plt.savefig(f"Plots/scatter_promoters_vs_enhancers_{col}_{cell}.pdf")
plt.close()






#--------------------------------------------------------------
# Archived
#--------------------------------------------------------------

def plot_observed_vs_expected(df,file):
    df_copy=df.copy()
    df_copy["num_rare_variants"]=np.log10(df_copy["num_rare_variants"]+1)
    df_copy["num_expected_variants"]=np.log10(df_copy["num_expected_variants"]+1)   
    # small dots
    sns.scatterplot(x="num_rare_variants", y="num_expected_variants", data=df_copy, hue="gc_content", s=6)
    # abline
    min_val=min(df_copy["num_rare_variants"].min(),df_copy["num_expected_variants"].min())
    max_val=max(df_copy["num_rare_variants"].max(),df_copy["num_expected_variants"].max())
    plt.plot([min_val,max_val],[min_val,max_val], color="black")
    plt.title(f'Observed v.s. Expected ({file})')
    plt.savefig(f"scatter_observed_expected_{file}.pdf")
    plt.close()


for dataset in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
    df_sub=df[df["dataset"]==dataset]
    sns.scatterplot(x="z", y="ism_cage", data=df_sub, hue="cooperativity",s=8)
    plt.title(f'{dataset}')
    plt.savefig(f"scatter_z_ism_cage_{dataset}.pdf")
    plt.close()
    sns.scatterplot(x="z", y="ism_starr", data=df_sub, hue="cooperativity",s=8)
    plt.title(f'{dataset}')
    plt.savefig(f"scatter_z_ism_starr_{dataset}.pdf")
    plt.close()



def plot_variable_vs_phylop_hue_by_gc(df,col,file):
    sns.scatterplot(x=col, y="241way_max", data=df, hue="gc_content")
    # calculate pearson correlation
    corr,pval = pearsonr(df[col], df["241way_max"])
    plt.text(0.3, 0.9, f"corr={corr.round(3)}\np={pval.round(3)}", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
    # annotate the name of the protein for all proteins in df_redundant, small text
    plt.title(f'{col} v.s. 241way_max ({file})')
    plt.savefig(f'scatter_{col}_hue_by_gc_241way_max_{file}.pdf')
    plt.close()


plot_variable_vs_phylop_hue_by_gc(df_sub,"z",file)






for re_type in ["promoters","enhancers"]:
    for coop in ["codependent","redundant"]:
        df_sub=df[df["re"]==re_type].reset_index(drop=True)
        df_sub=df_sub[df_sub["cooperativity"]==coop].reset_index(drop=True)
        corr_cage=pearsonr(df_sub["z"], df_sub["ism_cage"])
        corr_starr=pearsonr(df_sub["z"], df_sub["ism_starr"])
        logger.info(f"Correlation for {re_type} {coop}: cage {corr_cage}, starr {corr_starr}")




for cell in ["hepg2","k562"]:
    for coop in ["unknown"]:
        logger.info(f"Processing {cell} {coop}")
        df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
        df_sub=df_sub[df_sub["cooperativity"]==coop].reset_index(drop=True)
        ustat, pval = mannwhitneyu(df_sub[df_sub["re"]=="promoters"]["have_nonrare_variant"], df_sub[df_sub["re"]=="enhancers"]["have_nonrare_variant"])
        logger.info(f"U-statistic: {ustat}, p-value: {pval}")
        logger.info(f"Mean have_nonrare_variant: {df_sub[df_sub['re']=='promoters']['have_nonrare_variant'].mean()} {df_sub[df_sub['re']=='enhancers']['have_nonrare_variant'].mean()}")
        g = sns.jointplot(data=df_sub, x="have_nonrare_variant", y="241way_min", hue="re", marginal_kws={'common_norm': False})
        plt.subplots_adjust(top=0.9)  # Adjust the top to make space for the title
        g.fig.suptitle(f'{coop} tfs({cell})', fontsize=16)
        plt.savefig(f"Plots/joint_241way_min_have_nonrare_variant_promoters_vs_enhancers_{coop}_{cell}.pdf")
        plt.close()

