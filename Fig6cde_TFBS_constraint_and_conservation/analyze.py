import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


import matplotlib.cm as cm
import matplotlib.colors as mcolors




from scipy.stats import mannwhitneyu,pearsonr
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from seq_ops import SeqExtractor
from mutational_constraints import calc_constraint
from tf_cooperativity import assign_cooperativity 
from plotting import format_text




seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


#-------------------
# Helper functions
#-------------------




def add_tf_cooperativity(df,cell_line):
    if cell_line=="hepg2":
        thresh_redun=0.48
        thresh_codep=0.78
    elif cell_line=="k562":
        thresh_redun=0.43
        thresh_codep=0.80
    else:
        raise ValueError("cell line not found")
    # add tf ci
    df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_dhs.csv")
    df_coop=assign_cooperativity(df_coop,5,0.95,thresh_redun,thresh_codep)
    df_coop.rename(columns={"protein2":"protein"},inplace=True)
    df=pd.merge(df,df_coop,on="protein",how="inner")
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
    # TODO: change "have_common_variants" to "pct_common_variants"
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
# turn dataset into category, order: proximal_hepg2, distal_hepg2, proximal_k562, distal_k562
df_gnocchi["dataset"]=pd.Categorical(df_gnocchi["dataset"], categories=["proximal_hepg2","proximal_k562","distal_hepg2","distal_k562"], ordered=True)
df_gnocchi["dataset"] = df_gnocchi["dataset"].apply(format_text)




# calculate gnocchi constraint, and p value of 241way_max/mean/min for each cell
df_gnocchi_hepg2=aggregate_per_cell("hepg2")
df_gnocchi_k562=aggregate_per_cell("k562")





#----------------------------------
# 1. box plot "241way_max" by TF type
#----------------------------------


# Define color mapping
color_mapping = {
    "Independent": "white",
    "Redundant": "#1f77b4",  # cool color
    "Intermediate": "#7f7f7f",  # gray
    "Synergistic": "#d62728"  # warm color
}

for col in ["241way_max"]:
    plt.figure(figsize=(2.8, 2))
    # Thin frame
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    # Flier properties
    flierprops = {
        'marker': 'o',
        'markersize': 1,
        'markerfacecolor': 'gray',
        'markeredgecolor': 'gray',
        'alpha': 0.3,
    }
    # Custom palette mapping
    palette = {key: color_mapping[key] for key in df_gnocchi["cooperativity"].unique()}
    #
    sns.boxplot(data=df_gnocchi, x="dataset", y=col, hue="cooperativity", 
                palette=palette, linewidth=0.5, flierprops=flierprops)
    plt.xticks(fontsize=5, rotation=30)
    plt.yticks(fontsize=5)
    plt.xlabel("")
    plt.ylabel("max(PhyloP)", fontsize=7)
    # Anchor legend outside on the right
    plt.legend(
        title="TF type",
        fontsize=5,
        title_fontsize=5,
        loc="upper right",
    )
    plt.tight_layout()
    plt.savefig(f"box_{col}_vs_cooperativity.pdf")
    plt.close()





# calculate P values
df_hepg2=df_gnocchi[df_gnocchi["cell_line"]=="hepg2"].reset_index(drop=True)
df_hepg2_redundant=df_hepg2[df_hepg2["cooperativity"]=="Redundant"].reset_index(drop=True)
mannwhitneyu(df_hepg2_redundant[df_hepg2_redundant["re"]=="proximal"]["241way_max"], df_hepg2_redundant[df_hepg2_redundant["re"]=="distal"]["241way_max"])

df_hepg2_codependent=df_hepg2[df_hepg2["cooperativity"]=="Synergistic"].reset_index(drop=True)
mannwhitneyu(df_hepg2_codependent[df_hepg2_codependent["re"]=="proximal"]["241way_max"], df_hepg2_codependent[df_hepg2_codependent["re"]=="distal"]["241way_max"])

df_k562=df_gnocchi[df_gnocchi["cell_line"]=="k562"].reset_index(drop=True)
df_k562_redundant=df_k562[df_k562["cooperativity"]=="Redundant"].reset_index(drop=True)
mannwhitneyu(df_k562_redundant[df_k562_redundant["re"]=="proximal"]["241way_max"], df_k562_redundant[df_k562_redundant["re"]=="distal"]["241way_max"])

df_k562_codependent=df_k562[df_k562["cooperativity"]=="Sygernistic"].reset_index(drop=True)
mannwhitneyu(df_k562_codependent[df_k562_codependent["re"]=="proximal"]["241way_max"], df_k562_codependent[df_k562_codependent["re"]=="distal"]["241way_max"])

df_hepg2_proximal=df_hepg2[df_hepg2["re"]=="proximal"].reset_index(drop=True)
mannwhitneyu(df_hepg2_proximal[df_hepg2_proximal["cooperativity"]=="Redundant"]["241way_max"], df_hepg2_proximal[df_hepg2_proximal["cooperativity"]=="Synergistic"]["241way_max"])

df_hepg2_distal=df_hepg2[df_hepg2["re"]=="distal"].reset_index(drop=True)
mannwhitneyu(df_hepg2_distal[df_hepg2_distal["cooperativity"]=="Redundant"]["241way_max"], df_hepg2_distal[df_hepg2_distal["cooperativity"]=="Synergistic"]["241way_max"])

df_k562_proximal=df_k562[df_k562["re"]=="proximal"].reset_index(drop=True)
mannwhitneyu(df_k562_proximal[df_k562_proximal["cooperativity"]=="Redundant"]["241way_max"], df_k562_proximal[df_k562_proximal["cooperativity"]=="Synergistic"]["241way_max"])

df_k562_distal=df_k562[df_k562["re"]=="distal"].reset_index(drop=True)
mannwhitneyu(df_k562_distal[df_k562_distal["cooperativity"]=="Redundant"]["241way_max"], df_k562_distal[df_k562_distal["cooperativity"]=="Synergistic"]["241way_max"])






#----------------------------------
# 3. 241way_max vs cooperativity_index
#----------------------------------

# Get the list of unique datasets (files)
files = df_gnocchi["dataset"].unique()
x_var = "241way_max"

# Create a 2×2 grid of subplots with shared x and y axes
fig, axes = plt.subplots(2, 2, figsize=(3, 3), sharex=True, sharey=True)

for ax, file in zip(axes.flatten(), files):
    df_sub = df_gnocchi[df_gnocchi["dataset"] == file].reset_index(drop=True)
    df_sub = df_sub.dropna(subset=["cooperativity_index"]).reset_index(drop=True)
    # Scatter plot
    ax.scatter(df_sub[x_var], df_sub["cooperativity_index"], s=0.6,color="black", alpha=0.5)
    # Highlight CTCF if it exists
    df_ctcf = df_sub[df_sub["protein"].str.upper() == "CTCF"]
    if not df_ctcf.empty:
        ctcf_x = df_ctcf[x_var].values[0]
        ctcf_y = df_ctcf["cooperativity_index"].values[0]
        ax.plot(ctcf_x, ctcf_y, marker='o', color='red', markersize=3)
        ax.text(ctcf_x, ctcf_y, "CTCF", fontsize=5, color='red', va='bottom', ha='right')
    # Pearson correlation
    r, p_value = pearsonr(df_sub[x_var], df_sub["cooperativity_index"])
    # Frame and annotation
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.text(
        0.6, 0.8,
        f"r = {r:.2f}\np = {p_value:.2e}",
        fontsize=5,
        transform=ax.transAxes,
        verticalalignment='top',
        horizontalalignment='left',
    )
    ax.set_title(file, fontsize=7)
    ax.tick_params(axis='x', labelsize=5)
    ax.tick_params(axis='y', labelsize=5)


# Common axis labels
fig.text(0.5, 0.04, "max(PhyloP)", ha='center', fontsize=7)
fig.text(0.04, 0.5, "Synergy score", va='center', rotation='vertical', fontsize=7)

plt.tight_layout(rect=[0.05, 0.05, 1, 1])
plt.savefig(f"scatter_ci_vs_{x_var}.pdf", bbox_inches="tight")
plt.close()



#-------------------------------------------------------------------------------------------------------------------------------
# 4: separate by cooperativity, compare proximal and distal, scatter plot for those aggregate measurements (don't have p values)
#-------------------------------------------------------------------------------------------------------------------------------
# Define a custom color palette for cooperativity categories

custom_palette = {
    "Synergistic": "#d62728",  # warm color
    "Redundant": "#1f77b4",    # cool color
}



for cell in ["hepg2","k562"]:
    if cell=="hepg2":
        thresh_redun=0.48
        thresh_codep=0.78
    elif cell=="k562":
        thresh_redun=0.43
        thresh_codep=0.80
    else:
        raise ValueError("cell line not found")
    for col in ["z", "have_common_variants"]:
        # Subset the data to only constrained TFs for one cell line
        df_sub = df_gnocchi[df_gnocchi["cell_line"] == cell].reset_index(drop=True)
        df_pivot = df_sub.pivot(index="protein", columns="re", values=col).dropna().reset_index()
        if col == "have_common_variants":
            df_pivot["proximal"] *= 100
            df_pivot["distal"] *= 100
        #
        # Merge with stats
        if cell == "hepg2":
            df_pivot = df_pivot.merge(
                df_gnocchi_hepg2[["protein", "z", "p_val_241way_max", "p_val_241way_mean", "p_val_241way_min"]],
                on="protein", how="inner"
            )
        elif cell == "k562":
            df_pivot = df_pivot.merge(
                df_gnocchi_k562[["protein", "z", "p_val_241way_max", "p_val_241way_mean", "p_val_241way_min"]],
                on="protein", how="inner"
            )
        else:
            raise ValueError("Cell line not found")
        #
        # Keep only positive z-scores
        df_pivot = df_pivot[df_pivot["z"] > 0].reset_index(drop=True)
        #
        # Add cooperativity info
        df_pivot = add_tf_cooperativity(df_pivot, cell)
        df_pivot["cooperativity"] = df_pivot["cooperativity"].astype(str)
        #
        # Report stats
        df_redundant = df_pivot[df_pivot["cooperativity"] == "Redundant"].reset_index(drop=True)
        df_codependent = df_pivot[df_pivot["cooperativity"] == "Synergistic"].reset_index(drop=True)
        # "have_common_variants" margin: 0.01
        # "z" margin: 0.01
        logger.info(f"{cell.upper()} | {col} | Redundant TFs enhancers > promoters: {np.mean(df_redundant['distal'] > df_redundant['proximal']+0.1):.2f}")
        logger.info(f"{cell.upper()} | {col} | Synergistic TFs enhancers > promoters: {np.mean(df_codependent['distal'] > df_codependent['proximal']+0.1):.2f}")
        #
        # Special condition for outliers
        if col == "241way_max":
            df_outlier = df_pivot[(df_pivot["proximal"] > 4) & (df_pivot["distal"] > 4)].reset_index(drop=True)
            logger.info(f"Outlier proteins (enhancer>4 and promoter>4): {df_outlier.protein.tolist()}")
        #
        # === PLOT ===
        plt.figure(figsize=(2.6,2.2))
        ax = plt.gca()
        #
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)
        #
        # Redundant and Synergistic
        sns.scatterplot(
            x="proximal", y="distal",
            data=df_pivot[df_pivot["cooperativity"].isin(["Redundant", "Synergistic"])],
            hue="cooperativity",
            palette=custom_palette,
            s=5, legend=False
        )
        #
        # Intermediate TFs with color gradient
        intermediate_df = df_pivot[df_pivot["cooperativity"] == "Intermediate"]
        norm = mcolors.Normalize(vmin=thresh_redun, vmax=thresh_codep)
        cmap = cm.get_cmap("coolwarm")
        sc = ax.scatter(
            intermediate_df["proximal"], intermediate_df["distal"],
            c=intermediate_df["cooperativity_index"],
            cmap=cmap, norm=norm,
            marker="x", s=5, linewidths=0.5
        )
        #
        # Independent TFs
        sns.scatterplot(
            x="proximal", y="distal",
            data=df_pivot[df_pivot["cooperativity"] == "Independent"],s=9, label="Independent", marker="o", facecolors='none',edgecolor="black"
        )
        #
        # Diagonal reference line
        min_val = min(df_pivot["proximal"].min(), df_pivot["distal"].min())
        max_val = max(df_pivot["proximal"].max(), df_pivot["distal"].max())
        plt.plot([min_val, max_val], [min_val, max_val], color="black", linewidth=0.3)
        #
        # Custom legend handles
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', linestyle='None', label='Redundant',markerfacecolor=custom_palette["Redundant"], markersize=20),
            Line2D([0], [0], marker='o', color='w', linestyle='None', label='Synergistic',markerfacecolor=custom_palette["Synergistic"], markersize=20),
            Line2D([0], [0], marker='x', color='black', linestyle='None', label='Intermediate',markersize=20),
            Line2D([0], [0], marker='o', markerfacecolor='none', color='black', linestyle='None', label='Independent',markersize=20)
        ]
        ax.legend(handles=legend_elements, title="TF Type", fontsize=5, title_fontsize=5, markerscale=0.2)
        #
        # Add colorbar for cooperativity index
        cbar = plt.colorbar(sc, ax=ax)
        cbar.ax.tick_params(labelsize=5)
        cbar.set_label("Synergy score\n(Intermediate TFs)", fontsize=6)
        #
        # Labels and title
        if col == "z":
            plt.title("gnomAD constraint", fontsize=7)
            plt.xlabel("Promoter proximal", fontsize=7)
            plt.ylabel("Promoter distal", fontsize=7)
            plt.xticks(fontsize=5)
            plt.yticks(fontsize=5)
        elif col == "have_common_variants":
            plt.title("%TFBS with common SNP", fontsize=7)
            plt.xlabel("Promoter proximal (%)", fontsize=7)
            plt.ylabel("Promoter distal (%)", fontsize=7)
            plt.xticks(fontsize=5)
            plt.yticks(fontsize=5)
        else:
            plt.title(f"{col}", fontsize=7)
            plt.xlabel("Promoter proximal", fontsize=7)
            plt.ylabel("Promoter distal", fontsize=7)
            plt.xticks(fontsize=5)
            plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"proximal_vs_distal_{col}_{cell}.pdf", bbox_inches="tight")
        plt.close()