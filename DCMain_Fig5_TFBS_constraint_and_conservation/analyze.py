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



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


#-------------------
# Helper functions
#-------------------




def add_tf_cooperativity(df,cell_line):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt", header=None).iloc[:,0].tolist()
    df_tf_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}.csv")
    df_tf_coop.rename(columns={"protein2":"protein"}, inplace=True)
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    df=df.merge(df_tf_coop, on="protein", how="inner")
    return df




def read_file(file_suffix,cell_line,df=None):
    # if file_suffix is a string, read it
    if df is None:
        df=pd.read_csv(f'motif_info_thresh_500_dhs_{file_suffix}.csv')
    logger.info(df.shape)
    # remove first column
    # df.drop(df.columns[0], axis=1, inplace=True)
    track_nums=get_track_num(cell_line)
    track_nums_remove=[i for i in range(16) if i not in track_nums]
    cols_remove=[f"isa_track{i}" for i in track_nums_remove]
    df.drop(cols_remove, axis=1, inplace=True)
    mapper={0:"cage",1:"cage",2:"dhs",3:"dhs",4:"starr",5:"starr",6:"sure",7:"sure"}
    df.rename(columns={f"isa_track{i}":f"isa_{mapper[i]}" for i in track_nums}, inplace=True)
    df["max_af"] = df["gnomad_af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["min_af"] = df["gnomad_af"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["num_variants"] = df["gnomad_af"].apply(lambda x: len(str(x).split(":")))
    # does the loci have common variants?
    df["have_common_variants"] = df["gnomad_af"].apply(lambda x: any([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["gnomad_af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["241way_max"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_mean"] = df["phylop_241way"].apply(lambda x: np.mean([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["phylop_241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    df=add_tf_cooperativity(df,cell_line)
    return df


def summarize_by_tf(df):
    dataset=df["dataset"].iloc[0]
    df=df.groupby("protein").agg({"241way_max":"median",
                                  "241way_min":"median",
                                  "241way_mean":"median",
                                  "num_variants":"mean",
                                  "have_common_variants":"mean"
                                }).reset_index()
    df["dataset"]=dataset
    return df



#--------------------------------------------------------------
# Read in data
#--------------------------------------------------------------

def preprocess(file_name,cell_line,df=None):
    # if file_name is a string, read it
    df=read_file(file_name,cell_line,df)
    df_gnocchi=calc_constraint(df,seq_extractor)
    df_grouped=summarize_by_tf(df)
    df_gnocchi=df_gnocchi.merge(df_grouped, on="protein", how="inner")
    if "hepg2" in file_name:
        cell_line="hepg2"
    elif "k562" in file_name:
        cell_line="k562"
    else:
        raise ValueError("Invalid cell line")
    df_gnocchi=add_tf_cooperativity(df_gnocchi,cell_line)
    return df_gnocchi



# merge by protein
df_proximal_hepg2_gnocchi=preprocess("proximal_hepg2", "hepg2")
df_proximal_k562_gnocchi=preprocess("proximal_k562", "k562")

df_distal_hepg2_1=pd.read_csv("motif_info_thresh_500_dhs_distal_hepg2_v1.csv",nrows=1993215)
df_distal_hepg2_2=pd.read_csv("motif_info_thresh_500_dhs_distal_hepg2_v1.csv",skiprows=1993216,header=None)
df_distal_hepg2_2.drop(df_distal_hepg2_2.columns[29], axis=1, inplace=True)
df_distal_hepg2_2.columns=df_distal_hepg2_1.columns
df_distal_hepg2=pd.concat([df_distal_hepg2_1,df_distal_hepg2_2], axis=0, ignore_index=True)
df_distal_hepg2_gnocchi=preprocess("distal_hepg2","hepg2",df_distal_hepg2)


df_distal_k562_1=pd.read_csv("motif_info_thresh_500_dhs_distal_k562_v1.csv",nrows=938528)
df_distal_k562_2=pd.read_csv("motif_info_thresh_500_dhs_distal_k562_v1.csv",skiprows=938529,header=None)
df_distal_k562_2.drop(df_distal_k562_2.columns[29], axis=1, inplace=True)
df_distal_k562_2.columns=df_distal_k562_1.columns
df_distal_k562=pd.concat([df_distal_k562_1,df_distal_k562_2], axis=0, ignore_index=True)
df_distal_k562_gnocchi=preprocess("distal_k562","k562",df_distal_k562)




dict_files={"proximal_hepg2":df_proximal_hepg2_gnocchi,
            "proximal_k562":df_proximal_k562_gnocchi,
            "distal_hepg2":df_distal_hepg2_gnocchi,
            "distal_k562":df_distal_k562_gnocchi}




#-------------------------------------------------------------
# 1. scatter plot "241way_max" with "z" to see concordance, color by cooperativity
#-------------------------------------------------------------


for file in dict_files:
    # Calculate Pearson correlation and p-value
    r, p_value = pearsonr(dict_files[file]["z"], dict_files[file]["241way_max"])
    #
    plt.figure(figsize=(2.3, 2))
    # Thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    scatter = sns.scatterplot(
        x="z",
        y="241way_max",
        data=dict_files[file],
        hue="cooperativity_index",
        palette="coolwarm",
        s=8,
        legend=False
    )
    # Create a colorbar
    norm = plt.Normalize(dict_files[file]["cooperativity_index"].min(), dict_files[file]["cooperativity_index"].max())
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label("Cooperativity Index", fontsize=5)
    cbar.ax.tick_params(labelsize=5)
    #
    # Add Pearson r and p-value to the plot (top-left)
    plt.text(0.05, 0.95, f"Pearson R = {r:.2f}\np = {p_value:.2e}", fontsize=6, 
             transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left')
    #
    plt.tight_layout()
    plt.xlabel("z", fontsize=7)
    plt.ylabel("241way_max", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.savefig(f"Plots/scatter_z_241way_max_{file}.pdf")
    plt.close()



#----------------------------------
# 2. violin plot "241way_max" vs cooperativity
# remove unknown cooperativity
#----------------------------------

df=pd.concat([df_proximal_hepg2_gnocchi,df_proximal_k562_gnocchi,df_distal_hepg2_gnocchi,df_distal_k562_gnocchi], axis=0, ignore_index=True)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])
df=df[df["cooperativity"]!="unknown"].reset_index(drop=True)
# turn cooperativity into category, order: redudant, codependent
df["cooperativity"]=pd.Categorical(df["cooperativity"], categories=["redundant","codependent"], ordered=True)

# turn dataset into category, order: proximal_hepg2, distal_hepg2, proximal_k562, distal_k562
df["dataset"]=pd.Categorical(df["dataset"], categories=["proximal_hepg2","distal_hepg2","proximal_k562","distal_k562"], ordered=True)


plt.figure(figsize=(2.6,2.3))
# thin frame
plt.gca().spines['top'].set_linewidth(0.5)
plt.gca().spines['right'].set_linewidth(0.5)
plt.gca().spines['bottom'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)
sns.violinplot(data=df, 
               x="dataset", 
               y="241way_max", 
               cut=0,
               hue="cooperativity", 
               split=True,
               palette="coolwarm",
               inner="quartile",
               linewidth=0.5)
# x rotation 30 degrees
plt.xticks(rotation=30,fontsize=5)
plt.yticks(fontsize=5)
plt.xlabel("Regulatory element", fontsize=7)
plt.ylabel("241way_max", fontsize=7)
plt.legend(title="TF Cooperativity", fontsize=5, title_fontsize=5)
plt.tight_layout()
plt.savefig("violin_split_241way_max_vs_cooperativity.pdf")
plt.close()


# calculate P values
df_hepg2=df[df["cell_line"]=="hepg2"].reset_index(drop=True)
df_hepg2_redundant=df_hepg2[df_hepg2["cooperativity"]=="redundant"].reset_index(drop=True)
mannwhitneyu(df_hepg2_redundant[df_hepg2_redundant["re"]=="proximal"]["241way_max"], df_hepg2_redundant[df_hepg2_redundant["re"]=="distal"]["241way_max"])

df_hepg2_codependent=df_hepg2[df_hepg2["cooperativity"]=="codependent"].reset_index(drop=True)
mannwhitneyu(df_hepg2_codependent[df_hepg2_codependent["re"]=="proximal"]["241way_max"], df_hepg2_codependent[df_hepg2_codependent["re"]=="distal"]["241way_max"])

df_k562=df[df["cell_line"]=="k562"].reset_index(drop=True)
df_k562_redundant=df_k562[df_k562["cooperativity"]=="redundant"].reset_index(drop=True)
mannwhitneyu(df_k562_redundant[df_k562_redundant["re"]=="proximal"]["241way_max"], df_k562_redundant[df_k562_redundant["re"]=="distal"]["241way_max"])

df_k562_codependent=df_k562[df_k562["cooperativity"]=="codependent"].reset_index(drop=True)
mannwhitneyu(df_k562_codependent[df_k562_codependent["re"]=="proximal"]["241way_max"], df_k562_codependent[df_k562_codependent["re"]=="distal"]["241way_max"])

df_hepg2_proximal=df_hepg2[df_hepg2["re"]=="proximal"].reset_index(drop=True)
mannwhitneyu(df_hepg2_proximal[df_hepg2_proximal["cooperativity"]=="redundant"]["241way_max"], df_hepg2_proximal[df_hepg2_proximal["cooperativity"]=="codependent"]["241way_max"])

df_hepg2_distal=df_hepg2[df_hepg2["re"]=="distal"].reset_index(drop=True)
mannwhitneyu(df_hepg2_distal[df_hepg2_distal["cooperativity"]=="redundant"]["241way_max"], df_hepg2_distal[df_hepg2_distal["cooperativity"]=="codependent"]["241way_max"])

df_k562_proximal=df_k562[df_k562["re"]=="proximal"].reset_index(drop=True)
mannwhitneyu(df_k562_proximal[df_k562_proximal["cooperativity"]=="redundant"]["241way_max"], df_k562_proximal[df_k562_proximal["cooperativity"]=="codependent"]["241way_max"])

df_k562_distal=df_k562[df_k562["re"]=="distal"].reset_index(drop=True)
mannwhitneyu(df_k562_distal[df_k562_distal["cooperativity"]=="redundant"]["241way_max"], df_k562_distal[df_k562_distal["cooperativity"]=="codependent"]["241way_max"])






#----------------------------------
# 2. 241way_max vs cooperativity_index, and z vs cooperativity_index
# don't remove unknown cooperativity
#----------------------------------



df=pd.concat([df_proximal_hepg2_gnocchi,df_proximal_k562_gnocchi,df_distal_hepg2_gnocchi,df_distal_k562_gnocchi], axis=0, ignore_index=True)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])

variable_pairs = [
    ("z", "cooperativity_index"),
    ("241way_max", "cooperativity_index")
]

# Iterate over dict_files and variable pairs
for dataset_name, df in dict_files.items():
    for x_var, y_var in variable_pairs:
        # Calculate Pearson correlation
        r, p_value = pearsonr(df[x_var], df[y_var])
        #
        # Create scatter plot
        plt.figure(figsize=(2.3, 2))
        # thin frame
        plt.gca().spines['top'].set_linewidth(0.5)
        plt.gca().spines['right'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)
        plt.scatter(df[x_var], df[y_var], s=1)
        #
        # Add annotations for Pearson r and p-value
        plt.text(0.05, 0.95, f"r = {r:.2f}\np = {p_value:.2e}", fontsize=5, 
                 transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left')
        #
        # Add labels and title
        plt.xlabel(x_var, fontsize=7)
        plt.ylabel(y_var, fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        #
        # Save plot
        plt.savefig(f"Plots/scatter_{y_var}_vs_{x_var}_{dataset_name}.pdf")
        plt.close()








#--------------------------------------------------------------
# Analysis 3: separate by cooperativity, compare promoters and enhancers
#--------------------------------------------------------------


for cell in ["hepg2","k562"]:
    for col in ["z"]:   
        df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
        df_redundant=df_sub[df_sub["cooperativity"]=="redundant"].reset_index(drop=True)
        df_codependent=df_sub[df_sub["cooperativity"]=="codependent"].reset_index(drop=True)
        # pivot to compare z between promoters and enhancers
        df_redundant=df_redundant.pivot(index="protein", columns="re", values=col).dropna().reset_index()
        df_codependent=df_codependent.pivot(index="protein", columns="re", values=col).dropna().reset_index()
        # how many enhancers>promoters
        logger.info(f"Number of proteins where enhancers>promoters (redundant tfs): {np.sum(df_redundant['distal']>df_redundant['proximal'])}")
        logger.info(f"Number of proteins where enhancers>promoters (codependent tfs): {np.sum(df_codependent['distal']>df_codependent['proximal'])}")
        df_redundant["cooperativity"]="redundant"
        df_codependent["cooperativity"]="codependent"
        df_pivot=pd.concat([df_redundant,df_codependent], axis=0)
        #
        plt.figure(figsize=(2.3,2.2))
        # thin frame
        plt.gca().spines['top'].set_linewidth(0.5)
        plt.gca().spines['right'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)
        sns.scatterplot(x="proximal", y="distal", data=df_pivot, hue="cooperativity", s=5)
        # abline
        min_val=min(df_pivot["proximal"].min(),df_pivot["distal"].min())
        max_val=max(df_pivot["proximal"].max(),df_pivot["distal"].max())
        plt.plot([min_val,max_val],[min_val,max_val], color="black",linewidth=0.3)
        plt.legend(title="TF Cooperativity", fontsize=5, title_fontsize=5)
        plt.xlabel("Constraint at proximal", fontsize=7)
        plt.ylabel("Constraint at distal", fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"scatter_proximal_vs_distal_{col}_{cell}.pdf")
        plt.close()





