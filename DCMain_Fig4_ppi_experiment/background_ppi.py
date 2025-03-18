import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utils_ppi import read_pooled_baits_tf_pairs, mannwhitneyu_with_nan


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity




ci_suffix="pe"
redundancy_threshold=0.3 # 0.46
codependent_threshold=0.7 # 0.83



for mode in ["linearity_index","cooperativity_index"]:
    df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_{ci_suffix}.csv")
    df_coop=assign_cooperativity(df_coop,1,0.9,redundancy_threshold,codependent_threshold)
    df_coop=df_coop[df_coop["protein2"].isin(["BACH1","MAFG","IKZF1","RREB1","RFX5"])].reset_index(drop=True)

    # subset for proper background TFs (proteomics detectable), because TFs detected in proteomics tend to have lower ci
    proteins_background=pd.read_csv('Pd1_PPI_experiment/Nusinow_CCLE_K562proteomics.txt', sep='\t')["Gene_Symbol"].tolist()
    df_coop["in_proteomics"]=df_coop["protein1"].isin(proteins_background) & df_coop["protein2"].isin(proteins_background)
    # histogram of cooperativity index, split by "in_proteomics"
    sns.kdeplot(df_coop[df_coop["in_proteomics"]][mode],color="tab:blue",cut=0,label="In proteomics")
    sns.kdeplot(df_coop[~df_coop["in_proteomics"]][mode],color="tab:orange",cut=0,label="Not in proteomics")
    # add p value
    stat,p=mannwhitneyu_with_nan(df_coop[df_coop["in_proteomics"]][mode],df_coop[~df_coop["in_proteomics"]][mode])
    plt.text(0.5,0.5,f"p={p:.2e}",transform=plt.gca().transAxes)
    plt.xlabel(mode)
    plt.ylabel("Density")
    plt.title(f"Distribution of TF pair level {mode}")
    plt.legend()
    plt.savefig(f"background_5_baits_tf_pair_{mode}_{ci_suffix}.pdf")
    plt.close()



    # tf level
    df_tf=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_k562_{ci_suffix}.csv")
    df_tf["in_proteomics"]=df_tf["protein2"].isin(proteins_background)
    # histogram of cooperativity index, split by "in_proteomics"
    sns.kdeplot(df_tf[df_tf["in_proteomics"]][mode],color="tab:blue",cut=0,label="In proteomics")
    sns.kdeplot(df_tf[~df_tf["in_proteomics"]][mode],color="tab:orange",cut=0,label="Not in proteomics")
    # add p value
    stat,p=mannwhitneyu_with_nan(df_tf[df_tf["in_proteomics"]][mode],df_tf[~df_tf["in_proteomics"]][mode])
    plt.text(0.5,0.5,f"p={p:.2e}",transform=plt.gca().transAxes)
    plt.xlabel(mode)
    plt.ylabel("Density")
    plt.title(f"Distribution of TF level {mode}")
    plt.legend()
    plt.savefig(f"background_all_tf_{mode}_{ci_suffix}.pdf")
    plt.close()