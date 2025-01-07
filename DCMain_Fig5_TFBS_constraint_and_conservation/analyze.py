import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu,pearsonr
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from seq_ops import SeqExtractor
from mutational_constraints import calc_constraint

seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")

#-------------------
# Helper functions
#-------------------


def add_tf_cooperativity(df,cell_line):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df




def read_file(file_suffix):
    df=pd.read_csv(f'motif_info_thresh_500_dhs_{file_suffix}.csv')
    # remove first column
    df.drop(df.columns[0], axis=1, inplace=True)
    if 'hepg2' in file_suffix:
        cell_line="hepg2"
    elif 'k562' in file_suffix:
        cell_line="k562"
    track_nums=get_track_num(cell_line)
    track_nums_remove=[i for i in range(16) if i not in track_nums]
    cols_remove=[f"isa_track{i}" for i in track_nums_remove]
    df.drop(cols_remove, axis=1, inplace=True)
    mapper={0:"cage",1:"cage",2:"dhs",3:"dhs",4:"starr",5:"starr",6:"sure",7:"sure"}
    df.rename(columns={f"isa_track{i}":f"isa_{mapper[i]}" for i in track_nums}, inplace=True)
    df["max_af"] = df["gnomad_af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["min_af"] = df["gnomad_af"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["num_variants"] = df["gnomad_af"].apply(lambda x: len(str(x).split(":")))
    df["num_common_variants"] = df["gnomad_af"].apply(lambda x: sum([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["gnomad_af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["241way_max"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["phylop_241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    df=add_tf_cooperativity(df,cell_line)
    return df


def group_by_tf(df):
    dataset=df["dataset"].iloc[0]
    df=df.groupby("protein").agg({"241way_max":"median",
                                  "241way_min":"median",
                                  "num_variants":"sum",
                                  "num_common_variants":"sum"
                                }).reset_index()
    df["dataset"]=dataset
    return df



#--------------------------------------------------------------
# Read in data
#--------------------------------------------------------------
df_proximal_hepg2=read_file("proximal_hepg2")
df_proximal_k562=read_file("proximal_k562")
df_distal_hepg2=read_file("distal_hepg2")
df_distal_k562=read_file("distal_k562")

df_proximal_hepg2_gnocchi=calc_constraint(df_proximal_hepg2,seq_extractor)
df_proximal_k562_gnocchi=calc_constraint(df_proximal_k562,seq_extractor)
df_distal_hepg2_gnocchi=calc_constraint(df_distal_hepg2,seq_extractor)
df_distal_k562_gnocchi=calc_constraint(df_distal_k562,seq_extractor)


df_proximal_hepg2_grouped=group_by_tf(df_proximal_hepg2)
df_proximal_k562_grouped=group_by_tf(df_proximal_k562)
df_distal_hepg2_grouped=group_by_tf(df_distal_hepg2)
df_distal_k562_grouped=group_by_tf(df_distal_k562)



# merge by protein
df_proximal_hepg2_gnocchi=df_proximal_hepg2_gnocchi.merge(df_proximal_hepg2_grouped, on="protein", how="inner")
df_proximal_k562_gnocchi=df_proximal_k562_gnocchi.merge(df_proximal_k562_grouped, on="protein", how="inner")
df_distal_hepg2_gnocchi=df_distal_hepg2_gnocchi.merge(df_distal_hepg2_grouped, on="protein", how="inner")
df_distal_k562_gnocchi=df_distal_k562_gnocchi.merge(df_distal_k562_grouped, on="protein", how="inner")



# add tf cooperativity
df_proximal_hepg2_gnocchi=add_tf_cooperativity(df_proximal_hepg2_gnocchi,"hepg2")
df_proximal_k562_gnocchi=add_tf_cooperativity(df_proximal_k562_gnocchi,"k562")
df_distal_hepg2_gnocchi=add_tf_cooperativity(df_distal_hepg2_gnocchi,"hepg2")
df_distal_k562_gnocchi=add_tf_cooperativity(df_distal_k562_gnocchi,"k562")



# plot z vs 241way_max
sns.jointplot(x="z", y="241way_max", data=df_proximal_hepg2_gnocchi, hue="cooperativity", marginal_kws={'common_norm': False})
plt.savefig("joint_z_vs_241way_max_proximal_hepg2.pdf")
plt.close()


sns.jointplot(x="z", y="241way_max", data=df_proximal_k562_gnocchi, hue="cooperativity", marginal_kws={'common_norm': False})
plt.savefig("joint_z_vs_241way_max_proximal_k562.pdf")
plt.close()

sns.jointplot(x="z", y="241way_max", data=df_distal_hepg2_gnocchi, hue="cooperativity", marginal_kws={'common_norm': False})
plt.savefig("joint_z_vs_241way_max_distal_hepg2.pdf")
plt.close()

sns.jointplot(x="z", y="241way_max", data=df_distal_k562_gnocchi, hue="cooperativity", marginal_kws={'common_norm': False})
plt.savefig("joint_z_vs_241way_max_distal_k562.pdf")
plt.close()



# TODO: absolute number of variants



#--------------------------------------------------------------
# Analysis 2: separate by cooperativity, compare promoters and enhancers
#--------------------------------------------------------------


df=pd.concat([df_proximal_hepg2_gnocchi,df_proximal_k562_gnocchi,df_distal_hepg2_gnocchi,df_distal_k562_gnocchi], axis=0, ignore_index=True)
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
        ustat, pval = mannwhitneyu(df_sub[df_sub["re"]=="proximal"][xcol], df_sub[df_sub["re"]=="distal"][xcol])
        logger.info(f"U-statistic: {ustat}, p-value: {pval}")
        logger.info(f"Mean z score: {df_sub[df_sub['re']=='proximal'][xcol].mean()} {df_sub[df_sub['re']=='distal'][xcol].mean()}")
        g = sns.jointplot(data=df_sub, x=xcol, y=ycol,hue="re", marginal_kws={'common_norm': False})
        plt.subplots_adjust(top=0.9)  # Adjust the top to make space for the title
        g.fig.suptitle(f'{coop} tfs({cell})', fontsize=16)
        plt.savefig(f"Plots/joint_{ycol}_{xcol}_proximal_vs_distal_{coop}_{cell}.pdf")
        plt.close()






col="241way_max"
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
sns.scatterplot(x="proximal", y="distal", data=df_pivot, hue="cooperativity")
plt.title(f'{col}, {cell}')
# abline
min_val=min(df_pivot["proximal"].min(),df_pivot["distal"].min())
max_val=max(df_pivot["proximal"].max(),df_pivot["distal"].max())
plt.plot([min_val,max_val],[min_val,max_val], color="black")
plt.savefig(f"Plots/scatter_proximal_vs_distal_{col}_{cell}.pdf")
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
    sns.scatterplot(x="z", y="isa_cage", data=df_sub, hue="cooperativity",s=8)
    plt.title(f'{dataset}')
    plt.savefig(f"scatter_z_isa_cage_{dataset}.pdf")
    plt.close()
    sns.scatterplot(x="z", y="isa_starr", data=df_sub, hue="cooperativity",s=8)
    plt.title(f'{dataset}')
    plt.savefig(f"scatter_z_isa_starr_{dataset}.pdf")
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






for re_type in ["proximal","distal"]:
    for coop in ["codependent","redundant"]:
        df_sub=df[df["re"]==re_type].reset_index(drop=True)
        df_sub=df_sub[df_sub["cooperativity"]==coop].reset_index(drop=True)
        corr_cage=pearsonr(df_sub["z"], df_sub["isa_cage"])
        corr_starr=pearsonr(df_sub["z"], df_sub["isa_starr"])
        logger.info(f"Correlation for {re_type} {coop}: cage {corr_cage}, starr {corr_starr}")




for cell in ["hepg2","k562"]:
    for coop in ["unknown"]:
        logger.info(f"Processing {cell} {coop}")
        df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
        df_sub=df_sub[df_sub["cooperativity"]==coop].reset_index(drop=True)
        ustat, pval = mannwhitneyu(df_sub[df_sub["re"]=="proximal"]["have_nonrare_variant"], df_sub[df_sub["re"]=="distal"]["have_nonrare_variant"])
        logger.info(f"U-statistic: {ustat}, p-value: {pval}")
        logger.info(f"Mean have_nonrare_variant: {df_sub[df_sub['re']=='promoters']['have_nonrare_variant'].mean()} {df_sub[df_sub['re']=='enhancers']['have_nonrare_variant'].mean()}")
        g = sns.jointplot(data=df_sub, x="have_nonrare_variant", y="241way_min", hue="re", marginal_kws={'common_norm': False})
        plt.subplots_adjust(top=0.9)  # Adjust the top to make space for the title
        g.fig.suptitle(f'{coop} tfs({cell})', fontsize=16)
        plt.savefig(f"Plots/joint_241way_min_have_nonrare_variant_promoters_vs_enhancers_{coop}_{cell}.pdf")
        plt.close()

