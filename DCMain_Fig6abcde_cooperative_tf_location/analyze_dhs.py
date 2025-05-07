import pandas as pd
import numpy as np
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity 


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42




# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def read_df_add_tf_coop(prefix,file_name):
    df=pd.read_csv(f"{prefix}_{file_name}.csv",index_col=0)
    if "hepg2" in file_name:
        cell_line="hepg2"
        thresh_redun=0.48
        thresh_codep=0.78
    elif "k562" in file_name:
        cell_line="k562"
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



def count_tf(df,file_name):
    # count number of redundant and codependent TFs in each region
    df=df.groupby(["region","cooperativity"]).size().unstack(fill_value=0).reset_index()
    df.rename(columns={"Synergistic":"synergistic_tf_count",
                       "Redundant":"redundant_tf_count",
                       "Intermediate":"intermediate_tf_count",
                       "Independent":"independent_tf_count"},inplace=True)
    df["region_type"]=file_name
    return df



def assign_region_type(df):
    df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])
    df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[1])
    conditions = [
        (df['re'] == 'proximal') & (df['tissue_invariance'] == 'yes'),
        (df['re'] == 'proximal') & (df['tissue_invariance'] == 'no'),
        (df['re'] == 'distal') & (df['tissue_invariance'] == 'yes'),
        (df['re'] == 'distal') & (df['tissue_invariance'] == 'no'),
    ]
    choices = ['proximal_ti', 'proximal_ts', 'distal_ti', 'distal_ts']
    df["region_type"] = np.select(conditions, choices, default='distal_ts')
    df["region_type"] = pd.Categorical(
        df["region_type"],
        categories=["distal_ts", "distal_ti", "proximal_ts", "proximal_ti"],
        ordered=True
    )
    return df


def add_seq_info(df,file_name):
    df_seq_info=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/dhs_{file_name}.tsv",sep=" ")
    df_seq_info["region"]=df_seq_info["seqnames"]+":"+df_seq_info["start"].astype(str)+"-"+df_seq_info["end"].astype(str)
    # merge by region
    df=pd.merge(df,df_seq_info,on="region",how="inner")
    return df


prefix="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs"


#---------------------------------------
# DHS: boxplot for distribution of TF count
#---------------------------------------

def preprocess(prefix,file_name):
    df=read_df_add_tf_coop(prefix,file_name)
    df=count_tf(df,file_name)
    # add sequence info
    df=add_seq_info(df,file_name)
    df["dataset"]=file_name
    return df




df_proximal_hepg2=preprocess(prefix,"proximal_hepg2")
df_distal_hepg2=preprocess(prefix,"distal_hepg2")
df_proximal_k562=preprocess(prefix,"proximal_k562")
df_distal_k562=preprocess(prefix,"distal_k562")

df=pd.concat([df_proximal_hepg2,df_distal_hepg2,df_proximal_k562,df_distal_k562],axis=0)

df=assign_region_type(df)



# Neighbor pairs for Mann-Whitney U tests
neighbor_pairs = [
    ("distal_ts", "distal_ti"),
    ("distal_ti", "proximal_ts"),
    ("proximal_ts", "proximal_ti")
]

color_map = {
    "redundant_tf_count": "#1f77b4",
    "synergistic_tf_count": "#d62728",
    "intermediate_tf_count": "grey",
    "independent_tf_count": "black"
}



for cell_line in ["k562", "hepg2"]:
    df_sub = df[df["cell_line"] == cell_line].reset_index(drop=True)
    for col in ["redundant_tf_count", "synergistic_tf_count", "intermediate_tf_count", "independent_tf_count"]:
        col_color = color_map[col]
        logger.info(f"Plotting {col} for {cell_line}")
        plt.figure(figsize=(1.4, 2))
        ax = plt.gca()
        # thin frame, no spine on top and right  
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_linewidth(0.5)
        #
        # box plot with colored boxes and thin lines  
        if col!="linear_tf_count":
            sns.boxplot(
                x="region_type", 
                y=col, 
                data=df_sub, 
                color=col_color,  # this sets the default color
                boxprops={'facecolor': col_color, 'alpha': 0.3},  # add some transparency if desired
                whiskerprops={'color': "black"},
                capprops={'color': "black"},
                medianprops={'color': col_color},
                showfliers=False,
                linewidth=0.5
            )
        else:
            sns.boxplot(
                x="region_type", 
                y=col,
                data=df, 
                boxprops={'facecolor': 'none'},  
                whiskerprops={'color': "black"},
                capprops={'color': "black"},
                medianprops={'color': col_color},
                showfliers=True,  # Ensure outliers are shown
                linewidth=0.5,
                flierprops={'marker': 'o', 'markersize': 3, 'markerfacecolor': 'black', 'markeredgewidth': 0}  # Small black dots
            )
        #
        # Calculate and annotate Mann-Whitney U test p-values
        for pair in neighbor_pairs:
            group1 = df_sub[df_sub["region_type"] == pair[0]][col]
            group2 = df_sub[df_sub["region_type"] == pair[1]][col]
            #
            stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
            # Add p-value annotation at 95 percentile
            y_max = max(group1.quantile(0.95), group2.quantile(0.95))
            y_position = y_max * 1.05  # adjust multiplier as needed
            #
            x1 = df_sub["region_type"].cat.categories.get_loc(pair[0])
            x2 = df_sub["region_type"].cat.categories.get_loc(pair[1])
            #
            plt.plot([x1, x2], [y_position, y_position], lw=0.2, color='black')
            plt.text((x1 + x2) / 2, y_position * 1.05, f"p={p_value:.1e}",
                     ha='center', va='bottom', fontsize=5)
        #
        plt.xlabel("Region type", fontsize=7)
        plt.ylabel(f"{col}", fontsize=7)
        plt.xticks(fontsize=5,rotation=30)
        plt.yticks(fontsize=5)
        plt.tight_layout()
        plt.savefig(f"{col}_{cell_line}.pdf")
        plt.close()
























#---------------------------------------
# DHS: box plot for cooperativity index distribution of intermediate TFs (or all TFs)
#---------------------------------------

def preprocess2(prefix,file_name):
    df=read_df_add_tf_coop(prefix,file_name)
    df["dataset"]=file_name
    # add sequence info
    df=add_seq_info(df,file_name)
    # remove na
    df=df[~df["cooperativity_index"].isna()].reset_index(drop=True)
    df=assign_region_type(df)
    return df



df_proximal_hepg2=preprocess2(prefix,"proximal_hepg2")
df_distal_hepg2=preprocess2(prefix,"distal_hepg2")
df_proximal_k562=preprocess2(prefix,"proximal_k562")
df_distal_k562=preprocess2(prefix,"distal_k562")




df=pd.concat([df_proximal_hepg2,df_distal_hepg2,df_proximal_k562,df_distal_k562],axis=0)

# box plot cooperativity index by region type


def plot_synergy_by_region(df, cell, filename_prefix, neighbor_pairs):
    df_sub = df[df["cell_line"] == cell].reset_index(drop=True)
    plt.figure(figsize=(1.4, 2))
    ax = plt.gca()
    
    # Customize plot spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_linewidth(0.5)

    # Ensure region_type is categorical with correct order
    region_order = ["distal_ts", "distal_ti", "proximal_ts", "proximal_ti"]
    df_sub["region_type"] = pd.Categorical(df_sub["region_type"], categories=region_order, ordered=True)

    # Create boxplot
    sns.boxplot(
        x="region_type",
        y="cooperativity_index",
        data=df_sub,
        order=region_order,
        linewidth=0.5,
        boxprops={'facecolor': 'none'},
        whiskerprops={'color': 'black'},
        capprops={'color': 'black'}
    )
    plt.xticks(rotation=30)

    # Annotate Mann-Whitney U test p-values
    for pair in neighbor_pairs:
        group1 = df_sub[df_sub["region_type"] == pair[0]]["cooperativity_index"]
        group2 = df_sub[df_sub["region_type"] == pair[1]]["cooperativity_index"]
        stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
        y_max = max(group1.quantile(0.95), group2.quantile(0.95))
        y_position = y_max * 1.05
        x1 = region_order.index(pair[0])
        x2 = region_order.index(pair[1])
        plt.plot([x1, x2], [y_position, y_position], lw=0.2, color='black')
        plt.text((x1 + x2) / 2, y_position * 1.05, f"p={p_value:.1e}",
                 ha='center', va='bottom', fontsize=5)

    plt.xlabel("Region type", fontsize=7)
    plt.ylabel("Synergy score", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f"{filename_prefix}_{cell}.pdf")
    plt.close()









neighbor_pairs = [
    ("distal_ts", "distal_ti"),
    ("distal_ti", "proximal_ts"),
    ("proximal_ts", "proximal_ti")
]


# Plot for all TFs
for cell in ["hepg2", "k562"]:
    plot_synergy_by_region(df, cell, "tf_ci_by_region_type", neighbor_pairs)

# Filter and plot for Intermediate TFs
df_intermediate = df[df["cooperativity"] == "Intermediate"].reset_index(drop=True)
for cell in ["hepg2", "k562"]:
    plot_synergy_by_region(df_intermediate, cell, "intermediate_tf_ci_by_region_type", neighbor_pairs)





# nohup python3 analyze_dhs.py > analyze_dhs.log 2>&1 &















