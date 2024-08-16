import pandas as pd
import numpy as np
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer

df=pd.read_csv("tf_individual_effect_by_cell_type.csv")

def get_protein_list(df, ism):
       ism = f"dstat_ism_{ism}"
       df_sub = df[["protein","dataset", ism]]
       df_sub = df_sub.pivot(index='protein', columns='dataset', values=ism).reset_index()
       # does any row contain nan? No
       logger.info("Any row contain nan? {}".format(np.any(df_sub.isnull().any(axis=1))))
       df_sub["hepg2_k562"] = df_sub["hepg2"] * df_sub["k562"]
       df_sub["hepg2_k562_diff"] = df_sub["hepg2"] - df_sub["k562"]
       df_sub["hepg2_k562_sum"] = df_sub["hepg2"] + df_sub["k562"]
       proteins_hepg2_highlight=df_sub[(df_sub["hepg2_k562"] < 0) & (df_sub["hepg2_k562_diff"] > 0.3)]["protein"].tolist()
       proteins_k562_highlight = df_sub[(df_sub["hepg2_k562"] < 0) & (df_sub["hepg2_k562_diff"] < -0.3)]["protein"].tolist()
       proteins_both_highlight = df_sub[(df_sub["hepg2"] > 0.5 ) & (df_sub["k562"] > 0.5 ) ]["protein"].tolist()
       return proteins_hepg2_highlight, proteins_k562_highlight, proteins_both_highlight


def dict2df(this_dict):
       unique_tfs = sorted(list(set(sum(this_dict.values(), []))))
       df = pd.DataFrame(index=unique_tfs, columns=this_dict.keys())
       for col in this_dict:
              df[col] = df.index.isin(this_dict[col]).astype(int)
       return df


dict_hepg2_highlight = {"cage": [], "dhs": [], "starr": [], "sure": []}
dict_k562_highlight = {"cage": [], "dhs": [], "starr": [], "sure": []}
dict_both_highlight = {"cage": [], "dhs": [], "starr": [], "sure": []}


for ism in ["cage","dhs","starr","sure"]:
       proteins_hepg2_highlight, proteins_k562_highlight, proteins_both_highlight = get_protein_list(df, ism)
       dict_hepg2_highlight[ism] = proteins_hepg2_highlight
       dict_k562_highlight[ism] = proteins_k562_highlight
       dict_both_highlight[ism] = proteins_both_highlight


df_hepg2_highlight = dict2df(dict_hepg2_highlight)
df_k562_highlight = dict2df(dict_k562_highlight)
df_both_highlight = dict2df(dict_both_highlight)

df_hepg2_highlight.to_csv("tfs_hepg2_highlight.csv")
df_k562_highlight.to_csv("tfs_k562_highlight.csv")
df_both_highlight.to_csv("tfs_both_cells_highlight.csv")

# for df_hepg2_highlight
df_hepg2_highlight.sum(axis=1).value_counts() # 12
df_hepg2_highlight.shape # 76


df_k562_highlight.sum(axis=1).value_counts() # 29
df_k562_highlight.shape # 51


df_both_highlight.sum(axis=1).value_counts() # 18
df_both_highlight.shape # 54









#-----------
# Plotting
#-----------
def plot_dict(this_dict,title,output_file):
       df=dict2df(this_dict)
       plt.figure(figsize=(3, 16))
       my_plot = sns.heatmap(df, cmap="vlag", linewidths=0.5, annot=False)
       my_plot .set_xticklabels(my_plot .get_xticklabels(), rotation=90)
       plt.title(title)
       plt.xticks(rotation=90)
       plt.tight_layout()
       plt.savefig(output_file, dpi=300, bbox_inches='tight')
       plt.close()


plot_dict(dict_hepg2_highlight,"TFs highlighted in HepG2","Plots/tfs_hepg2_highlight.pdf")
plot_dict(dict_k562_highlight, "TFs highlighted in K562","Plots/tfs_k562_highlight.pdf")
plot_dict(dict_both_highlight, "TFs highlighted in both HepG2 and K562","Plots/tfs_both_highlight.pdf")



#--------------------------------------
# Add information of RNA expression
#--------------------------------------

def get_expr(tfs, df_expr):
       logger.info(f"Number of TFs: {len(tfs)}")
       df_expr_sub=df_expr[df_expr["HGNC.symbol"].isin(tfs)]
       logger.info(f"Number of TFs found in expression data: {df_expr_sub.shape[0]}")
       df_expr_sub=df_expr_sub.set_index("HGNC.symbol")
       df_expr_sub=df_expr_sub.loc[tfs]
       tfs_expr=np.log(df_expr_sub["TPM"].values+1)
       return tfs_expr

df_expr_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Curate_motif_annot/TF_expression_hepg2.csv")
df_expr_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Curate_motif_annot/TF_expression_k562.csv")


def get_tf_expr_table(df):
       tfs = df.index.tolist()
       tfs=split_dimer(tfs)
       tfs_expr_hepg2=get_expr(tfs,df_expr_hepg2)
       tfs_expr_k562=get_expr(tfs,df_expr_k562)
       df_expr=pd.DataFrame({"TF":tfs,"HepG2":tfs_expr_hepg2,"K562":tfs_expr_k562})
       # set TF as index
       df_expr=df_expr.set_index("TF")
       return df_expr


df_hepg2_highlight_expr=get_tf_expr_table(df_hepg2_highlight)
df_hepg2_highlight_expr["TF_type"]="Highlighted by HepG2"
df_k562_highlight_expr=get_tf_expr_table(df_k562_highlight)
df_k562_highlight_expr["TF_type"]="Highlighted by K562"
df_both_highlight_expr=get_tf_expr_table(df_both_highlight)
df_both_highlight_expr["TF_type"]="Highlighted by both"


# mann-whitney test
from scipy.stats import mannwhitneyu
mannwhitneyu(df_hepg2_highlight_expr["HepG2"],df_k562_highlight_expr["K562"])
mannwhitneyu(df_k562_highlight_expr["HepG2"],df_k562_highlight_expr["K562"])
mannwhitneyu(df_both_highlight_expr["HepG2"],df_both_highlight_expr["K562"])

mannwhitneyu(df_k562_highlight_expr["K562"],df_both_highlight_expr["K562"])
mannwhitneyu(df_hepg2_highlight_expr["HepG2"],df_both_highlight_expr["HepG2"])




df_plot=pd.concat([df_hepg2_highlight_expr,df_k562_highlight_expr,df_both_highlight_expr])

# change to long format for violin plot
df_plot=df_plot.reset_index()
df_plot=pd.melt(df_plot,id_vars=["TF","TF_type"],value_vars=["HepG2","K562"],var_name="cell_type",value_name="TPM")


plt.figure(figsize=(8, 6))
sns.boxplot(x="TF_type",y="TPM",hue="cell_type",data=df_plot)
plt.title("TF expression in HepG2 and K562")
plt.savefig("Plots/tf_expression.pdf",dpi=300,bbox_inches='tight')
plt.close()