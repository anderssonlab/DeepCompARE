import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity

#-------------------------------
# calculate cell type specificity
#-------------------------------


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42





threshold_dict={"hepg2":{"pe":[0.3,0.7],"dhs":[0.48,0.78]},
                "k562":{"pe":[0.3,0.7],"dhs":[0.44,0.81]}}


suffix="pe"

for cell_line in ["hepg2","k562"]:
    df_dispersion=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/gtex.dispersionEstimates.tab",sep="\t")
    #
    df_tf=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_{suffix}.csv")
    df_tf=assign_cooperativity(df_tf,
                            5,0.95,
                            threshold_dict[cell_line][suffix][0],
                            threshold_dict[cell_line][suffix][1])
    # remove independent
    df_tf=df_tf[df_tf["cooperativity"]!="Independent"].reset_index(drop=True)
    # turn df_tf["cooperativity"] to categorical
    df_tf["cooperativity"]=pd.Categorical(df_tf["cooperativity"],categories=["Redundant","Intermediate","Synergistic"],ordered=True)
    # merge with df_dispersion
    df_tf=df_tf.merge(df_dispersion,left_on="protein2",right_on="symbol",how="inner")
    df_tf.drop_duplicates(subset="protein2",inplace=True)
    #
    
    # scatter plot for gini and cooperativity index
    plt.figure(figsize=(2.3,2.3))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(y="gini",x="cooperativity_index",data=df_tf,hue="cooperativity",palette={"Intermediate":"gray","Synergistic":"#d62728" ,"Redundant":"#1f77b4"},s=5,alpha=0.5)
    # add pearson correlation and p
    r,p=pearsonr(df_tf["gini"],df_tf["cooperativity_index"])
    plt.text(0.1,1.2,f"Pearson R={r:.2f}\nP={p:.2e}",fontsize=5)
    # Add horizontal median lines for each cooperativity group
    group_colors = {"Intermediate":"gray", "Synergistic":"#d62728", "Redundant":"#1f77b4"}
    median_x_positions = {
        "Redundant": 0.15,
        "Intermediate": 0.5,
        "Synergistic": 0.85
    }

    for group, x_center in median_x_positions.items():
        median_val = df_tf[df_tf["cooperativity"] == group]["gini"].median()
        plt.hlines(y=median_val, xmin=x_center - 0.05, xmax=x_center + 0.05,
               color=group_colors[group], linewidth=2.0)
    plt.xlabel("Synergy score", fontsize=7)
    plt.ylabel("Cell type specificity (Gini index)", fontsize=7)
    # legend
    plt.legend(title="TF type",fontsize=5, title_fontsize=5, markerscale=0.5)
    # xticks should include 
    plt.xticks([0,0.3,0.7,1],fontsize=5)
    plt.yticks([0.2,0.4,0.6,0.8,1,1.2,1.4],fontsize=5)
    plt.title(cell_line, fontsize=7)
    plt.tight_layout()
    plt.savefig(f"gini_vs_ci_{cell_line}_{suffix}.pdf")
    plt.close()
