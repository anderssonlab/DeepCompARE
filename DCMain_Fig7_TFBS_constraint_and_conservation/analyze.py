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


# add mediator info to compare conservation


def add_tf_cooperativity(df,cell_line):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}_dhs.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}_dhs.txt", header=None).iloc[:,0].tolist()
    df_tf_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_dhs.csv")
    df_tf_coop=df_tf_coop[df_tf_coop["c_sum"]>1].reset_index(drop=True)
    df_tf_coop.rename(columns={"protein2":"protein"}, inplace=True)
    df["cooperativity"]="intermediate"
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
    # does the loci have common variants?
    df["have_common_variants"] = df["gnomad_af"].apply(lambda x: np.any([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["gnomad_af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["241way_max"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_mean"] = df["phylop_241way"].apply(lambda x: np.mean([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["phylop_241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=this_file
    df=add_tf_cooperativity(df,cell_line)
    return df




def calculate_pval_proximal_vs_distal(df_orig,col):
    # remove rows with nan in col
    df=df_orig.dropna(subset=[col]).reset_index(drop=True)
    dict_results={
        "protein":[],
        "p_val":[],
    }
    for protein in df.protein.unique():
        df_tf=df[df["protein"]==protein].reset_index(drop=True)
        df_proximal=df_tf[df_tf["re"]=="proximal"].reset_index(drop=True)
        df_distal=df_tf[df_tf["re"]=="distal"].reset_index(drop=True)
        if len(df_proximal)<10 or len(df_distal)<10:
            continue
        # mannwhitneyu test
        _, pval = mannwhitneyu(df_proximal[col], df_distal[col])
        dict_results["protein"].append(protein)
        dict_results["p_val"].append(pval)
    return pd.DataFrame(dict_results)





def summarize_by_tf(df):
    return df.groupby("protein").agg({"241way_max":"median",
                                      "241way_min":"median",
                                      "241way_mean":"median",
                                      "have_common_variants":"mean"}).reset_index()




def preprocess(this_file):
    df=read_file(this_file)
    if "hepg2" in this_file:
        cell_line="hepg2"
    elif "k562" in this_file:
        cell_line="k562"
    else:
        raise ValueError("cell line not found")
    df_gnocchi=calc_constraint(df,seq_extractor)
    df_grouped=summarize_by_tf(df)
    df_gnocchi=df_gnocchi.merge(df_grouped, on="protein", how="inner")
    df_gnocchi=add_tf_cooperativity(df_gnocchi,cell_line)
    df_gnocchi["dataset"]=this_file
    return df_gnocchi



def aggregate_per_cell(cell_line):
    df_proximal=read_file(f"proximal_{cell_line}")
    df_distal=read_file(f"distal_{cell_line}")
    df_cell=pd.concat([df_proximal,df_distal], axis=0, ignore_index=True)
    df_cell["re"]=df_cell["dataset"].apply(lambda x: x.split("_")[0])
    df_241way_max_pval=calculate_pval_proximal_vs_distal(df_cell,"241way_max")
    df_241way_mean_pval=calculate_pval_proximal_vs_distal(df_cell,"241way_mean")
    df_241way_min_pval=calculate_pval_proximal_vs_distal(df_cell,"241way_min")
    df_241way_max_pval.rename(columns={"p_val":"p_val_241way_max"}, inplace=True)
    df_241way_mean_pval.rename(columns={"p_val":"p_val_241way_mean"}, inplace=True)
    df_241way_min_pval.rename(columns={"p_val":"p_val_241way_min"}, inplace=True)
    df_res=df_241way_max_pval.merge(df_241way_mean_pval, on="protein", how="inner")
    df_res=df_res.merge(df_241way_min_pval, on="protein", how="inner")
    # add constraint
    df_constraint=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_{cell_line}_dhs.csv")
    df_res=df_res.merge(df_constraint, on="protein", how="inner")
    return df_res






#--------------------------------------------------------------
# Read in data
#--------------------------------------------------------------

# calculate gnocchi constraint for each file
df_proximal_hepg2_gnocchi=preprocess("proximal_hepg2")
df_proximal_k562_gnocchi=preprocess("proximal_k562")
df_distal_hepg2_gnocchi=preprocess("distal_hepg2")
df_distal_k562_gnocchi=preprocess("distal_k562")


df_gnocchi=pd.concat([df_proximal_hepg2_gnocchi,df_proximal_k562_gnocchi,df_distal_hepg2_gnocchi,df_distal_k562_gnocchi], axis=0, ignore_index=True)
df_gnocchi["cell_line"]=df_gnocchi["dataset"].apply(lambda x: x.split("_")[1])
df_gnocchi["re"]=df_gnocchi["dataset"].apply(lambda x: x.split("_")[0])
# turn cooperativity into category, order: redudant, codependent
df_gnocchi["cooperativity"]=pd.Categorical(df_gnocchi["cooperativity"], categories=["redundant","intermediate","codependent"], ordered=True)
# turn dataset into category, order: proximal_hepg2, distal_hepg2, proximal_k562, distal_k562
df_gnocchi["dataset"]=pd.Categorical(df_gnocchi["dataset"], categories=["proximal_hepg2","distal_hepg2","proximal_k562","distal_k562"], ordered=True)





# calculate gnocchi constraint, and p value of 241way_max/mean/min for each cell
df_gnocchi_hepg2=aggregate_per_cell("hepg2")
df_gnocchi_k562=aggregate_per_cell("k562")


#-------------------------------------------------------------
# 1. scatter plot "241way_max" with "z" to see concordance, color by cooperativity
#-------------------------------------------------------------


for file in df_gnocchi["dataset"].unique():
    df_sub=df_gnocchi[df_gnocchi["dataset"]==file].reset_index(drop=True)
    # Calculate Pearson correlation and p-value
    r, p_value = pearsonr(df_sub["z"], df_sub["241way_max"])
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
        data=df_sub,
        hue="cooperativity_index",
        palette="coolwarm",
        s=8,
        legend=False
    )
    # Create a colorbar
    norm = plt.Normalize(df_sub["cooperativity_index"].min(), df_sub["cooperativity_index"].max())
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
# 2. box plot "241way_max" vs cooperativity
#----------------------------------


for col in ["241way_max","241way_mean","241way_min"]:
    plt.figure(figsize=(2.6,2.3))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    flierprops = {
        'marker': 'o',        # marker style
        'markersize': 1,      # marker size
        'markerfacecolor': 'gray',  # marker face color
        'markeredgecolor': 'gray'  # marker edge color
    }
    sns.boxplot(data=df_gnocchi, x="dataset", y=col, hue="cooperativity", palette="coolwarm",linewidth=0.5, flierprops=flierprops)
    # x rotation 30 degrees
    plt.xticks(rotation=30,fontsize=5)
    plt.yticks(fontsize=5)
    plt.xlabel("Regulatory element", fontsize=7)
    plt.ylabel(col, fontsize=7)
    plt.legend(title="TF Cooperativity", fontsize=5, title_fontsize=5, loc="upper right")
    plt.tight_layout()
    plt.savefig(f"box_{col}_vs_cooperativity.pdf")
    plt.close()



# calculate P values
df_hepg2=df_gnocchi[df_gnocchi["cell_line"]=="hepg2"].reset_index(drop=True)
df_hepg2_redundant=df_hepg2[df_hepg2["cooperativity"]=="redundant"].reset_index(drop=True)
mannwhitneyu(df_hepg2_redundant[df_hepg2_redundant["re"]=="proximal"]["241way_max"], df_hepg2_redundant[df_hepg2_redundant["re"]=="distal"]["241way_max"])

df_hepg2_codependent=df_hepg2[df_hepg2["cooperativity"]=="codependent"].reset_index(drop=True)
mannwhitneyu(df_hepg2_codependent[df_hepg2_codependent["re"]=="proximal"]["241way_max"], df_hepg2_codependent[df_hepg2_codependent["re"]=="distal"]["241way_max"])

df_k562=df_gnocchi[df_gnocchi["cell_line"]=="k562"].reset_index(drop=True)
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
# don't remove intermediate cooperativity
#----------------------------------


# Iterate over dict_files and variable pairs
for file in df_gnocchi["dataset"].unique():
    df_sub = df_gnocchi[df_gnocchi["dataset"] == file].reset_index(drop=True)
    for x_var in["z","241way_max"]:
        # Calculate Pearson correlation
        r, p_value = pearsonr(df_sub[x_var], df_sub["cooperativity_index"])
        #
        # Create scatter plot
        plt.figure(figsize=(2.3, 2))
        # thin frame
        plt.gca().spines['top'].set_linewidth(0.5)
        plt.gca().spines['right'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)
        plt.scatter(df_sub[x_var], df_sub["cooperativity_index"], s=1)
        #
        # Add annotations for Pearson r and p-value
        plt.text(0.6, 0.25, f"r = {r:.2f}\np = {p_value:.2e}", fontsize=5, 
                 transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left')
        #
        # Add labels and title
        plt.xlabel(x_var, fontsize=7)
        plt.ylabel("Cooperativity index", fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        #
        # Save plot
        plt.savefig(f"scatter_ci_vs_{x_var}_{file}.pdf")
        plt.close()








#-------------------------------------------------------------------------------------------------------------------------------
# 4: separate by cooperativity, compare proximal and distal, scatter plot for those aggregate measurements (don't have p values)
#-------------------------------------------------------------------------------------------------------------------------------
# Define a custom color palette for cooperativity categories

from matplotlib.lines import Line2D

custom_palette = {
    "codependent": "orangered",  # warm color
    "redundant": "dodgerblue",   # cool color
}
for cell in ["hepg2", "k562"]:
    # for col in ["z", "have_common_variants"]:
    for col in ["z"]:
        # Subset the data to only constrained TFs for one cell line
        df_sub = df_gnocchi[df_gnocchi["cell_line"] == cell].reset_index(drop=True)
        df_pivot = df_sub.pivot(index="protein", columns="re", values=col).dropna().reset_index()
        # merge 
        if cell=="hepg2":
            df_pivot=df_pivot.merge(df_gnocchi_hepg2[["protein","z","p_val_241way_max","p_val_241way_mean","p_val_241way_min"]], on="protein", how="inner")
        elif cell=="k562":
            df_pivot=df_pivot.merge(df_gnocchi_k562[["protein","z","p_val_241way_max","p_val_241way_mean","p_val_241way_min"]], on="protein", how="inner")
        else:
            raise ValueError("cell line not found")
        #
        df_pivot=df_pivot[df_pivot["z"]>0].reset_index(drop=True)
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
        sns.scatterplot(x="proximal", y="distal", data=df_pivot[df_pivot["cooperativity"]!="intermediate"], hue="cooperativity", s=5, palette=custom_palette)
        sns.scatterplot(x="proximal", y="distal", data=df_pivot[df_pivot["cooperativity"]=="intermediate"], hue="cooperativity_index", s=5, palette="coolwarm",legend=False,marker="x")
        min_val = min(df_pivot["proximal"].min(), df_pivot["distal"].min())
        max_val = max(df_pivot["proximal"].max(), df_pivot["distal"].max())
        plt.plot([min_val, max_val], [min_val, max_val], color="black", linewidth=0.3)
        #
        dummy_handle = Line2D([], [],
                      marker='x',
                      color='black',
                      markersize=5,
                      linestyle='None',
                      markeredgewidth=0.5,
                      label='intermediate')
        #   
        # Get current handles and labels, and append the dummy handle:
        handles, labels = ax.get_legend_handles_labels()
        handles.append(dummy_handle)
        labels.append("intermediate")
        ax.legend(handles=handles, labels=labels, title="TF Cooperativity", fontsize=5, title_fontsize=5,markerscale=0.5)
        #
        plt.xlabel(f"{col} at proximal", fontsize=7)
        plt.ylabel(f"{col} at distal", fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"scatter_proximal_vs_distal_{col}_{cell}.pdf")
        plt.close()



#--------------------------------------------------------------
# 5: separate by cooperativity, compare proximal and distal
#--------------------------------------------------------------

for cell in ["hepg2", "k562"]:
    for col in ["241way_max", "241way_mean", "241way_min"]:
        # Subset the data to only constrained TFs for one cell line
        df_sub = df_gnocchi[df_gnocchi["cell_line"] == cell].reset_index(drop=True)
        df_pivot = df_sub.pivot(index="protein", columns="re", values=col).dropna().reset_index()
        # merge 
        if cell=="hepg2":
            df_pivot=df_pivot.merge(df_gnocchi_hepg2[["protein","z","p_val_241way_max","p_val_241way_mean","p_val_241way_min"]], on="protein", how="inner")
        elif cell=="k562":
            df_pivot=df_pivot.merge(df_gnocchi_k562[["protein","z","p_val_241way_max","p_val_241way_mean","p_val_241way_min"]], on="protein", how="inner")
        else:
            raise ValueError("cell line not found")
        #
        df_pivot=df_pivot[df_pivot["z"]>0].reset_index(drop=True)
        # Add cooperativity information
        df_pivot = add_tf_cooperativity(df_pivot, cell)
        # determine alpha: if f"p_val_{col}"<0.05, alpha=1, else alpha=0.3
        df_pivot["alpha"]=df_pivot.apply(lambda x: 1 if x[f"p_val_{col}"]<0.05 else 0.3, axis=1)
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
        sns.scatterplot(
            x="proximal", 
            y="distal", 
            data=df_pivot[df_pivot["cooperativity"] != "intermediate"],
            hue="cooperativity", 
            s=5, 
            palette=custom_palette, 
            alpha=df_pivot[df_pivot["cooperativity"] != "intermediate"]["alpha"]
        )
        #
        sns.scatterplot(
            x="proximal", 
            y="distal", 
            data=df_pivot[df_pivot["cooperativity"] == "intermediate"],
            hue="cooperativity_index", 
            s=5, 
            palette="coolwarm",
            legend=False,
            marker="x", 
            alpha=df_pivot[df_pivot["cooperativity"] == "intermediate"]["alpha"]
        )
        min_val = min(df_pivot["proximal"].min(), df_pivot["distal"].min())
        max_val = max(df_pivot["proximal"].max(), df_pivot["distal"].max())
        plt.plot([min_val, max_val], [min_val, max_val], color="black", linewidth=0.3)
        #
        dummy_handle = Line2D([], [],
                      marker='x',
                      color='black',
                      markersize=5,
                      linestyle='None',
                      markeredgewidth=0.5,
                      label='intermediate')
        #   
        # Get current handles and labels, and append the dummy handle:
        handles, labels = ax.get_legend_handles_labels()
        handles.append(dummy_handle)
        labels.append("intermediate")
        ax.legend(handles=handles, labels=labels, title="TF Cooperativity", fontsize=5, title_fontsize=5,markerscale=0.5)
        plt.xlabel(f"{col} at proximal", fontsize=7)
        plt.ylabel(f"{col} at distal", fontsize=7)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"scatter_proximal_vs_distal_{col}_{cell}.pdf")
        plt.close()

