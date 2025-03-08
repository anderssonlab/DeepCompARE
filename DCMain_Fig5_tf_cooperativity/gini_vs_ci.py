import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer
#-------------------------------
# calculate cell type specificity
#-------------------------------



suffix="pe"

for cell_line in ["hepg2","k562"]:
    df_dispersion=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/gtex.dispersionEstimates.tab",sep="\t")
    #
    df_tf=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_{suffix}.csv")
    df_tf=df_tf[df_tf["c_sum"]>1].reset_index(drop=True)
    # merge with df_dispersion
    df_tf=df_tf.merge(df_dispersion,left_on="protein2",right_on="symbol",how="inner")
    df_tf.drop_duplicates(subset="protein2",inplace=True)
    #
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}_{suffix}.txt",header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}_{suffix}.txt",header=None).iloc[:,0].tolist()
    df_tf["tf_type"]="intermediate"
    df_tf.loc[df_tf["protein2"].isin(tfs_codependent),"tf_type"]="codependent"
    df_tf.loc[df_tf["protein2"].isin(tfs_redundant),"tf_type"]="redundant"
    # scatter plot for gini and cooperativity index
    plt.figure(figsize=(2.5,2.5))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(y="gini",x="cooperativity_index",data=df_tf,hue="tf_type",palette={"intermediate":"gray","codependent":"orangered","redundant":"dodgerblue"},s=5)
    # add pearson correlation and p
    r,p=pearsonr(df_tf["gini"],df_tf["cooperativity_index"])
    plt.text(0.05,0.7,f"r={r:.2f}\np={p:.2e}",fontsize=5)
    plt.xlabel("Cooperativity index", fontsize=7)
    plt.ylabel("Cell type specificity (Gini)", fontsize=7)
    # legend
    plt.legend(title="TF type",fontsize=5, title_fontsize=5)
    # xticks should include 
    plt.xticks(fontsize=5)
    plt.yticks([0.25,0.5,0.75,1],fontsize=5)
    plt.tight_layout()
    plt.savefig(f"gini_vs_ci_{cell_line}_{suffix}.pdf")
    plt.close()
