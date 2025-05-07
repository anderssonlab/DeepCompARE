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

prefix="Pd2_E1E2P1P2_motif_info/motif_info_thresh_500"


def read_df_add_tf_coop(prefix,file_name):
    df=pd.read_csv(f"{prefix}_{file_name}_k562.csv",index_col=0)
    # add tf ci
    df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_k562_dhs.csv")
    df_coop=assign_cooperativity(df_coop,5,0.95,0.43,0.80)
    df_coop.rename(columns={"protein2":"protein"},inplace=True)
    df=pd.merge(df,df_coop,on="protein",how="inner")
    df["region_type"]=file_name
    return df



def count_tf(df,file_name):
    # count number of redundant and codependent TFs in each region
    df=df.groupby(["region","cooperativity"]).size().unstack(fill_value=0).reset_index()
    df.rename(columns={"Synergistic":"synergistic_tf_count",
                       "Redundant":"redundant_tf_count",
                       "Intermediate":"intermediate_tf_count",
                       "Independent":"independent_tf_count"
                       },inplace=True)
    df["region_type"]=file_name
    return df



#---------------------------------------
# E1E2P1P2: boxplot for distribution of TF count
#---------------------------------------

def preprocess(prefix,file_name):
    df=read_df_add_tf_coop(prefix,file_name)
    df=count_tf(df,file_name)
    return df




df_p1=preprocess(prefix,"P1")
df_p2=preprocess(prefix,"P2")
df_e1=preprocess(prefix,"E1")
df_e2=preprocess(prefix,"E2")





df=pd.concat([df_p1,df_p2,df_e1,df_e2],axis=0)
df["region_type"]=pd.Categorical(df["region_type"],categories=["E1","E2","P2","P1"],ordered=True)




# Neighbor pairs for Mann-Whitney U tests
neighbor_pairs = [
    ("E1", "E2"),
    ("E2", "P2"),
    ("P2", "P1")
]

color_map = {
    "redundant_tf_count": "#1f77b4",
    "synergistic_tf_count": "#d62728",
    "intermediate_tf_count": "grey",
    "independent_tf_count": "black"
}



for col in ["redundant_tf_count", "synergistic_tf_count", "intermediate_tf_count", "independent_tf_count"]:
    col_color = color_map[col]
    plt.figure(figsize=(1.6, 2.3))
    ax = plt.gca()
    # thin frame, no top and right
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    #
    # if col != "independent_tf_count":
    sns.boxplot(
        x="region_type", 
        y=col, 
        data=df, 
        color=col_color,  # this sets the default color
        boxprops={'facecolor': col_color, 'alpha': 0.3},  # add some transparency if desired
        whiskerprops={'color': "black"},
        capprops={'color': "black"},
        medianprops={'color': col_color},
        showfliers=False,
        linewidth=0.5,
    )
    # if col == "independent_tf_count":
    #     sns.boxplot(
    #         x="region_type", 
    #         y="independent_tf_count",
    #         data=df, 
    #         boxprops={'facecolor': 'none'},  
    #         whiskerprops={'color': "black"},
    #         capprops={'color': "black"},
    #         medianprops={'color': col_color},
    #         showfliers=True,  # Ensure outliers are shown
    #         linewidth=0.5,
    #         flierprops={'marker': 'o', 'markersize': 3, 'markerfacecolor': 'black', 'markeredgewidth': 0}  # Small black dots
    #     )

    # Calculate and annotate Mann-Whitney U test p-values
    for pair in neighbor_pairs:
        group1 = df[df["region_type"] == pair[0]][col]
        group2 = df[df["region_type"] == pair[1]][col]
        stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
        # Add p-value annotation at 95 percentile
        y_max = max(group1.quantile(0.95), group2.quantile(0.95))
        y_position = y_max * 1.05  # adjust multiplier as needed
        #
        x1 = df["region_type"].cat.categories.get_loc(pair[0])
        x2 = df["region_type"].cat.categories.get_loc(pair[1])
        #
        plt.plot([x1, x2], [y_position, y_position], lw=0.2, color='black')
        plt.text((x1 + x2) / 2, y_position * 1.05, f"p={p_value:.2e}",
                    ha='center', va='bottom', fontsize=5)
        
    plt.xlabel("Region type", fontsize=7)
    plt.ylabel(f"{col} count", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f"e1e2p1p2_{col}.pdf")
    plt.close()


















#---------------------------------------
# E1E2P1P2: box plot for cooperativity index distribution of intermediate TFs
#---------------------------------------
df_p1=read_df_add_tf_coop(prefix,"P1")
df_p2=read_df_add_tf_coop(prefix,"P2")
df_e1=read_df_add_tf_coop(prefix,"E1")
df_e2=read_df_add_tf_coop(prefix,"E2")

df=pd.concat([df_p1,df_p2,df_e1,df_e2],axis=0)

# select cooperativity=="intermediate"
df=df[df["cooperativity"]=="Intermediate"].reset_index(drop=True)
# df=df[df["cooperativity"]!="Independent"].reset_index(drop=True)

df["region_type"]=pd.Categorical(df["region_type"],categories=["E1","E2","P2","P1"],ordered=True)


# box plot cooperativity index by region type
neighbor_pairs = [
    ("E1", "E2"),
    ("E2", "P2"),
    ("P2", "P1")
]


plt.figure(figsize=(1.6, 2.3))
ax = plt.gca()
# thin frame, no top and right
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for spine in ax.spines.values():
    spine.set_linewidth(0.5)

sns.boxplot(
    x="region_type",
    y="cooperativity_index",
    data=df,
    linewidth=0.5,
    boxprops={'facecolor': 'none'},
    whiskerprops={'color': 'black'},
    capprops={'color': 'black'},
    flierprops={'marker': 'o', 'markersize': 1}
)
#
# Calculate and annotate Mann-Whitney U test p-values
for pair in neighbor_pairs:
    group1 = df[df["region_type"] == pair[0]]["cooperativity_index"]
    group2 = df[df["region_type"] == pair[1]]["cooperativity_index"]
    stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
    # Determine y position for annotation (using 95th percentile as reference)
    y_max = max(group1.quantile(0.95), group2.quantile(0.95))
    y_position = y_max * 1.05  # Adjust as needed
    # Get x-axis positions based on the categorical order
    x1 = df["region_type"].cat.categories.get_loc(pair[0])
    x2 = df["region_type"].cat.categories.get_loc(pair[1])
    # Draw a line connecting the two groups
    plt.plot([x1, x2], [y_position, y_position], lw=0.2, color='black')
    # Add p-value text; adjust multiplier to position text closer to the line if needed
    plt.text((x1 + x2) / 2, y_position * 1.05, f"p={p_value:.2e}",
                ha='center', va='bottom', fontsize=5)

plt.xlabel("Region type", fontsize=7)
plt.ylabel("Synergy score", fontsize=7)
plt.xticks(fontsize=5)
plt.yticks(fontsize=5)
plt.tight_layout()
plt.savefig(f"e1e2p1p2_intermediate_tf_ci_by_region_type.pdf")
plt.close()







# nohup python3 re_ci.py > re_ci.out &















