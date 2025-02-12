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
from scipy.stats import ks_2samp

seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


#-------------------
# Helper functions
#-------------------


def add_tf_cooperativity(df,cell_line):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}_dhs.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}_dhs.txt", header=None).iloc[:,0].tolist()
    df_tf_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_dhs.csv")
    df_tf_coop.rename(columns={"protein2":"protein"}, inplace=True)
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    df=df.merge(df_tf_coop, on="protein", how="inner")
    return df




def read_file(this_file):
    # if this_file is a string, read it
    prefix_motif_info="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs"
    df=pd.read_csv(f'{prefix_motif_info}_{this_file}.csv')
    if "hepg2" in this_file:
        cell_line="hepg2"
    elif "k562" in this_file:
        cell_line="k562"
    else:
        raise ValueError("cell line not found")
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
    df["have_common_variants"] = df["gnomad_af"].apply(lambda x: np.any([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["gnomad_af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["241way_max"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_mean"] = df["phylop_241way"].apply(lambda x: np.mean([float(i) for i in str(x).split(":")]))
    df["241way_median"] = df["phylop_241way"].apply(lambda x: np.median([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["phylop_241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=this_file
    df=add_tf_cooperativity(df,cell_line)
    return df


def summarize_by_tf(df):
    df=df.groupby("protein").agg({"241way_max":"median",
                                  "241way_min":"median",
                                  "241way_mean":"median",
                                  "241way_median":"median",
                                  "have_common_variants":"mean",
                                }).reset_index()
    return df


def preprocess(this_file,cell_line):
    df=read_file(this_file)
    df_gnocchi=calc_constraint(df,seq_extractor)
    df_grouped=summarize_by_tf(df)
    df_gnocchi=df_gnocchi.merge(df_grouped, on="protein", how="inner")
    df_gnocchi=add_tf_cooperativity(df_gnocchi,cell_line)
    df_gnocchi["dataset"]=this_file
    return df_gnocchi






#--------------------------------------------------------------
# Read in data
#--------------------------------------------------------------



# merge by protein
df_proximal_hepg2_gnocchi=preprocess("proximal_hepg2", "hepg2")
df_proximal_k562_gnocchi=preprocess("proximal_k562", "k562")
df_distal_hepg2_gnocchi=preprocess("distal_hepg2", "hepg2")
df_distal_k562_gnocchi=preprocess("distal_k562", "k562")




dict_files={"proximal_hepg2":df_proximal_hepg2_gnocchi,
            "proximal_k562":df_proximal_k562_gnocchi,
            "distal_hepg2":df_distal_hepg2_gnocchi,
            "distal_k562":df_distal_k562_gnocchi}




# get somewhat constrained TFs
dict_tfs_constrained={}

for cell in ["hepg2","k562"]:
    df_proximal=read_file(f"proximal_{cell}")
    df_distal=read_file(f"distal_{cell}")
    df_all=pd.concat([df_proximal,df_distal], axis=0, ignore_index=True)
    df_all_gnocchi=calc_constraint(df_all,seq_extractor)
    df_constrained=df_all_gnocchi[df_all_gnocchi["z"]>0].reset_index(drop=True)
    dict_tfs_constrained[cell]=df_constrained["protein"].unique().tolist()



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
    plt.savefig(f"scatter_z_241way_max_{file}.pdf")
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


col="241way_max"
plt.figure(figsize=(2.6,2.3))
# thin frame
plt.gca().spines['top'].set_linewidth(0.5)
plt.gca().spines['right'].set_linewidth(0.5)
plt.gca().spines['bottom'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)
sns.violinplot(data=df, 
               x="dataset", 
               y=col,
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
plt.ylabel(col, fontsize=7)
plt.legend(title="TF Cooperativity", fontsize=5, title_fontsize=5)
plt.tight_layout()
plt.savefig(f"violin_split_{col}_vs_cooperativity.pdf")
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
# 3. 241way_max vs cooperativity_index, and z vs cooperativity_index
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
        plt.text(0.6, 0.25, f"r = {r:.2f}\np = {p_value:.2e}", fontsize=5, 
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
        plt.savefig(f"scatter_{y_var}_vs_{x_var}_{dataset_name}.pdf")
        plt.close()








#-------------------------------------------------------------------------------------------------------------------------------
# 4: separate by cooperativity, compare promoters and enhancers, scatter plot for those aggregate measurements (don't have p values)
#-------------------------------------------------------------------------------------------------------------------------------
# Define a custom color palette for cooperativity categories
custom_palette = {
    "codependent": "orangered",  # warm color
    "redundant": "dodgerblue",   # cool color
    "unknown": "#C0C0C0"
}
for cell in ["hepg2", "k562"]:
    for col in ["z", "have_common_variants"]:
        # Subset the data to only constrained TFs for one cell line
        df_sub = df[df["cell_line"] == cell].reset_index(drop=True)
        df_sub = df_sub[df_sub["protein"].isin(dict_tfs_constrained[cell])].reset_index(drop=True)
        df_pivot = df_sub.pivot(index="protein", columns="re", values=col).dropna().reset_index()
        # Add cooperativity information
        df_pivot = add_tf_cooperativity(df_pivot, cell)
        df_redundant = df_pivot[df_pivot["cooperativity"] == "redundant"].reset_index(drop=True)
        df_codependent = df_pivot[df_pivot["cooperativity"] == "codependent"].reset_index(drop=True)
        logger.info(f"Number of proteins where enhancers>promoters (redundant tfs): {np.mean(df_redundant['distal'] > df_redundant['proximal'])}")
        logger.info(f"Number of proteins where enhancers>promoters (codependent tfs): {np.mean(df_codependent['distal'] > df_codependent['proximal'])}")
        # Create the plot
        plt.figure(figsize=(2.3, 2.2))
        ax = plt.gca()
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        #
        # Use the custom palette in the scatterplot
        sns.scatterplot(x="proximal", y="distal", data=df_pivot, hue="cooperativity", s=5, palette=custom_palette)
        #
        min_val = min(df_pivot["proximal"].min(), df_pivot["distal"].min())
        max_val = max(df_pivot["proximal"].max(), df_pivot["distal"].max())
        plt.plot([min_val, max_val], [min_val, max_val], color="black", linewidth=0.3)
        #
        plt.legend(title="TF Cooperativity", fontsize=5, title_fontsize=5)
        plt.xlabel(f"{col} at proximal", fontsize=7)
        plt.ylabel(f"{col} at distal", fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"scatter_proximal_vs_distal_{col}_{cell}.pdf")
        plt.close()



#--------------------------------------------------------------
# 5: separate by cooperativity, compare promoters and enhancers， volcano plot
#--------------------------------------------------------------

df_proximal_hepg2=read_file(f"proximal_hepg2")
df_distal_hepg2=read_file(f"distal_hepg2")
df_proximal_k562=read_file(f"proximal_k562")
df_distal_k562=read_file(f"distal_k562")

df=pd.concat([df_proximal_hepg2,df_distal_hepg2,df_proximal_k562,df_distal_k562], axis=0, ignore_index=True)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])



dict_results={
    "protein":[],
    "diff_241way_max":[],
    "p_241way_max":[],
    "diff_241way_min":[],
    "p_241way_min":[]
}




def calculate_volcano(cell):
    df_cell=df[df["cell_line"]==cell].reset_index(drop=True)
    for protein in dict_tfs_constrained[cell]:
        df_tf=df_cell[df_cell["protein"]==protein].reset_index(drop=True)
        df_proximal=df_tf[df_tf["re"]=="proximal"].reset_index(drop=True)
        df_distal=df_tf[df_tf["re"]=="distal"].reset_index(drop=True)
        # KS test
        _, p_max = ks_2samp(df_proximal["241way_max"], df_distal["241way_max"])
        diff_max=df_proximal["241way_max"].median()-df_distal["241way_max"].median()
        dict_results["protein"].append(protein)
        dict_results["diff_241way_max"].append(diff_max)
        dict_results["p_241way_max"].append(p_max)
        # KS test
        _, p_min = ks_2samp(df_proximal["241way_min"], df_distal["241way_min"])
        diff_min=df_proximal["241way_min"].median()-df_distal["241way_min"].median()
        dict_results["diff_241way_min"].append(diff_min)
        dict_results["p_241way_min"].append(p_min)
        print(f"{protein}, {diff_max}, {p_max}, {diff_min}, {p_min}")
    return pd.DataFrame(dict_results)




df_volcano_hepg2=calculate_volcano("hepg2")
# add cooperativity
df_volcano_hepg2=add_tf_cooperativity(df_volcano_hepg2,"hepg2")

df_volcano_k562=calculate_volcano("k562")
# add cooperativity
df_volcano_k562=add_tf_cooperativity(df_volcano_k562,"k562")




# plot volcano plot

df_plot=df_volcano_hepg2
col="241way_min"
# log10 transform p values
df_plot["-log10(p)"]=df_plot[f"p_{col}"].apply(lambda x: -np.log10(x))
# clip p values at 30
df_plot["-log10(p)"]=df_plot["-log10(p)"].clip(upper=30)
# scatter
plt.figure(figsize=(2.3, 2.2))
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(0.5)







