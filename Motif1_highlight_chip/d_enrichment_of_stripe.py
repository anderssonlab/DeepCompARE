import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2_contingency

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")

#------------------------
# plot for separate file
#------------------------
def plot_enrichment_separate_file(ks_file,cell_type,title):
    df=pd.read_csv(f"summary_and_ks_test_{ks_file}.csv")
    df["fdr"]=multipletests(df.feat_imp_p_val, method="fdr_bh")[1]
    df["significant"]=df["fdr"]<0.05
    df_stripe=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_potential_{cell_type}.csv") 
    df_stripe["is_stripe"]=df_stripe["colocolized_by_t2"]>df_stripe.shape[0]/2
    df=df.merge(df_stripe[["tf2","is_stripe"]], left_on="protein", right_on="tf2", how="inner")
    contingency_table=pd.crosstab(df["significant"], df["is_stripe"])
    chi2, p, dof, expected = chi2_contingency(contingency_table) 
    sns.scatterplot(data=df, x="mean_feat_imp", y="mean_feat_imp_true", hue="is_stripe")
    plt.text(0.95, 0.05, f"chi2 p={p:.3f}", ha='center', va='center', transform=plt.gca().transAxes)
    plt.xlabel("Mean feature importance")
    plt.ylabel("Mean feature importance (ChIP=True)")
    # add diagonal line, limit is max(max(x),max(y))
    lim=max(df["mean_feat_imp"].max(),df["mean_feat_imp_true"].max())
    plt.plot([0,lim],[0,lim],color="black",linestyle="--")
    plt.title(title)
    plt.savefig(f"Plots/mean_feat_imp_{ks_file}.pdf",dpi=300)
    plt.close()
    
# plot_enrichment_separate_file("enhancers_hepg2","HepG2","Enhancer HepG2")
# plot_enrichment_separate_file("promoters_hepg2","HepG2","Promoters HepG2")
# plot_enrichment_separate_file("enhancers_k562","K562","Enhancer K562")
# plot_enrichment_separate_file("promoters_k562","K562","Promoters K562")



def plot_enrichment_all_files():
    df1=pd.read_csv(f"summary_and_ks_test_enhancers_hepg2.csv")
    df2=pd.read_csv(f"summary_and_ks_test_promoters_hepg2.csv")
    df3=pd.read_csv(f"summary_and_ks_test_enhancers_k562.csv")
    df4=pd.read_csv(f"summary_and_ks_test_promoters_k562.csv")
    df=pd.concat([df1,df2,df3,df4])
    df["fdr"]=multipletests(df.feat_imp_p_val, method="fdr_bh")[1]
    # for rows with same "protein", select the row with the lowest fdr
    df=df.groupby("protein").agg({"mean_feat_imp":"mean","mean_feat_imp_true":"mean","fdr":"min"}).reset_index()
    df["significant"]=df["fdr"]<0.05
    
    df_stripe1=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_potential_HepG2.csv") 
    df_stripe2=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_potential_K562.csv")
    df_stripe1["is_stripe"]=df_stripe1["colocolized_by_t2"]>df_stripe1.shape[0]/2
    df_stripe2["is_stripe"]=df_stripe2["colocolized_by_t2"]>df_stripe2.shape[0]/2
    df_stripe=pd.concat([df_stripe1,df_stripe2])
    # for rows with same "tf2", "is_stripe" is the "OR" of all rows
    df_stripe=df_stripe.groupby("tf2").agg({"is_stripe":"max"}).reset_index() 

    df=df.merge(df_stripe[["tf2","is_stripe"]], left_on="protein", right_on="tf2", how="inner")
    contingency_table=pd.crosstab(df["significant"], df["is_stripe"])
    chi2, p, dof, expected = chi2_contingency(contingency_table) 
    sns.scatterplot(data=df, x="mean_feat_imp", y="mean_feat_imp_true", hue="is_stripe")
    plt.text(0.95, 0.05, f"chi2 p={p:.3f}", ha='center', va='center', transform=plt.gca().transAxes)
    plt.xlabel("Mean feature importance")
    plt.ylabel("Mean feature importance (ChIP=True)")
    # add diagonal line, limit is max(max(x),max(y))
    lim=max(df["mean_feat_imp"].max(),df["mean_feat_imp_true"].max())
    plt.plot([0,lim],[0,lim],color="black",linestyle="--")
    plt.title("Mean feature importance (all files)")
    plt.savefig(f"Plots/mean_feat_imp_all.pdf",dpi=300)
    plt.close()
    
plot_enrichment_all_files()
    