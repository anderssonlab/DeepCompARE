import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2_contingency

from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")


df_cage=pd.read_csv("Pd1_ks_tests/summary_ks_test_cage.csv")
df_cage["fdr"]=multipletests(df_cage.ism_motif_p_val, method="fdr_bh")[1]
df_cage=df_cage[df_cage.fdr<0.05].reset_index(drop=True)
df_cage["file_name"]="cage"
np.sum(df_cage.ism_motif_d_stat<0) # 19/245


df_dhs=pd.read_csv("Pd1_ks_tests/summary_ks_test_dhs.csv")
df_dhs["fdr"]=multipletests(df_dhs.ism_motif_p_val, method="fdr_bh")[1]
df_dhs=df_dhs[df_dhs.fdr<0.05].reset_index(drop=True)
df_dhs["file_name"]="dhs"
np.sum(df_dhs.ism_motif_d_stat<0) # 61/250

df_starr=pd.read_csv("Pd1_ks_tests/summary_ks_test_starr.csv")
df_starr["fdr"]=multipletests(df_starr.ism_motif_p_val, method="fdr_bh")[1]
df_starr=df_starr[df_starr.fdr<0.05].reset_index(drop=True)
df_starr["file_name"]="starr"
np.sum(df_starr.ism_motif_d_stat<0) # 118/232

df_sure=pd.read_csv("Pd1_ks_tests/summary_ks_test_sure.csv")
df_sure["fdr"]=multipletests(df_sure.ism_motif_p_val, method="fdr_bh")[1]
df_sure=df_sure[df_sure.fdr<0.05].reset_index(drop=True)
df_sure["file_name"]="sure"
np.sum(df_sure.ism_motif_d_stat<0) # 49/231


df=pd.concat([df_cage,df_dhs,df_starr,df_sure]).reset_index(drop=True)
df=df[df.motif_score_p_val>0.05]

df["fdr"]=multipletests(df.ism_motif_p_val, method="fdr_bh")[1]

np.sum(df.motif_score_p_val<0.05) # 56/660, ignore-able
np.sum(df.fdr<0.05) # 239/660 TF with fdr<0.05
np.sum(df.feat_imp_d_stat>0) # 508/660
pearsonr(df_cage.ism_motif_d_stat, df_cage.occupancy) # 0.19, 0.002
pearsonr(df_dhs.ism_motif_d_stat, df_dhs.occupancy) # 0.04, 0.45
pearsonr(df_starr.ism_motif_d_stat, df_starr.occupancy) # -0.08, 0.19
pearsonr(df_sure.ism_motif_d_stat, df_sure.occupancy) # 0.05, 0.40



#---------------------------------
# Analysis1: scatter plot, color dot by dataset
#---------------------------------

texts = []
categories = df['file_name'].unique()

# Generate a color palette with seaborn
palette = sns.color_palette("husl", len(categories))
palette_dict = dict(zip(categories, palette))

for category in categories:
    # Filter data for the current category
    category_data = df[df['file_name'] == category]
    # Plot each category with its corresponding color
    scatter = plt.scatter(category_data['feat_imp_d_stat'], category_data['conditional_occupancy'],
                          color=palette_dict[category], label=category,s=1)
    # Add text for each point in the current category
    for line in range(0, category_data.shape[0]):
        text = plt.text(category_data['feat_imp_d_stat'].iloc[line], category_data['conditional_occupancy'].iloc[line], 
                        category_data['protein'].iloc[line], horizontalalignment='left', 
                        size='xx-small', color=palette_dict[category])
        texts.append(text)

# Adjust texts to avoid overlap
adjust_text(texts)
# add correlation annotation at bottom right
plt.text(x=plt.xlim()[0] + (plt.xlim()[1] - plt.xlim()[0]) * 0.02, 
         y=plt.ylim()[1] - (plt.ylim()[1] - plt.ylim()[0]) * 0.1, 
         s=f"correlation={corr:.3f}", 
         fontsize=10, ha='right', va='bottom')
plt.legend()
plt.xlabel("KS d statistics of feat_imp")
plt.ylabel("Conditional occupancy")
plt.title("KS d statistics of feat_imp and conditional occupancy")
plt.savefig("Plots_ks_test/ks_dstat_featimp_and_occupancy_scatter.pdf",dpi=300)
plt.close()




#---------------------------------
# Analysis2: Robin's heatmap
#---------------------------------
# seaborn heatmap
# convert df to wide format, with file_name as rows, protein as columns, and feat_imp_d_stat as values
df_dstat=df.pivot(index="file_name", columns="protein", values="feat_imp_d_stat")
df_pval=df.pivot(index="file_name", columns="protein", values="ism_motif_p_val")
df_fdr=df.pivot(index="file_name", columns="protein", values="fdr")

df_dstat_filled = df_dstat.fillna(-1)
df_pval_filled = df_pval.fillna(-1)
df_fdr_filled = df_fdr.fillna(-1)

# Perform hierarchical clustering on the filled data
row_linkage = linkage(squareform(pdist(df_dstat_filled)), method='average')
col_linkage = linkage(squareform(pdist(df_dstat_filled.T)), method='average')

# Get the leaves order
df_dstat_filled = df_dstat_filled.iloc[leaves_list(row_linkage), :]
df_dstat_filled = df_dstat_filled.iloc[:, leaves_list(col_linkage)]
df_pval_filled = df_pval_filled.iloc[leaves_list(row_linkage), :]
df_pval_filled = df_pval_filled.iloc[:, leaves_list(col_linkage)]
df_fdr_filled = df_fdr_filled.iloc[leaves_list(row_linkage), :]
df_fdr_filled = df_fdr_filled.iloc[:, leaves_list(col_linkage)]

df_dstat_for_plotting = df_dstat_filled.replace(-1, np.nan)
df_pval_for_plotting = df_pval_filled.replace(-1, np.nan)
df_fdr_for_plotting = df_fdr_filled.replace(-1, np.nan)

plt.figure(figsize=(50, 3))
sns.heatmap(df_dstat_for_plotting, cmap='coolwarm', cbar_kws={'label': 'KS d statistics'},xticklabels=True)
for i in range(df_fdr_for_plotting.shape[0]):  # Rows
    for j in range(df_fdr_for_plotting.shape[1]):  # Columns
        if df_fdr_for_plotting.iloc[i, j] < 0.05:
            plt.text(j+0.5, i+0.5, '*', ha='center', va='center', color='red')
plt.subplots_adjust(bottom=0.5)
plt.title("KS d statistics of feat_imp")
plt.savefig("Plots_ks_test/ks_dstat_featimp_heatmap_fdr.pdf",dpi=300)
plt.close()

plt.figure(figsize=(50, 3))
sns.heatmap(df_dstat_for_plotting, cmap='coolwarm', cbar_kws={'label': 'KS d statistics'},xticklabels=True)
for i in range(df_pval_for_plotting.shape[0]):  # Rows
    for j in range(df_pval_for_plotting.shape[1]):  # Columns
        if df_pval_for_plotting.iloc[i, j] < 0.05:
            plt.text(j+0.5, i+0.5, '*', ha='center', va='center', color='red')
plt.subplots_adjust(bottom=0.5)
plt.title("KS d statistics of feat_imp")
plt.savefig("Plots_ks_test/ks_dstat_featimp_heatmap_pval.pdf",dpi=300)
plt.close()

# Fun stuff
# All significant: TCF7, GABPA, ZNF384, TCF7L2
# Promoter biased: TBP, NR3C1
# K562 biased: ELK1::SREBF2

#------------------------------------------------------------------
# Analysis3: correspondance between ism_motif_p_val and remove_context_p_val
#------------------------------------------------------------------

pearsonr(df.feat_imp_d_stat, df.remove_context_d_stat) # 0.17, 3e-5
    
# get contingincy table for the two p-values
df["ism_motif_p_val_bin"]=pd.cut(df.ism_motif_p_val, bins=[0,0.05,1], labels=["significant","non-significant"])
df["remove_context_p_val_bin"]=pd.cut(df.remove_context_p_val, bins=[0,0.05,1], labels=["significant","non-significant"])
contingency_table=pd.crosstab(df.ism_motif_p_val_bin, df.remove_context_p_val_bin)

# chi2 test

chi2, p, dof, expected = chi2_contingency(contingency_table) # p=0.001

    