import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect.csv")
df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["RE"]=df["dataset"].apply(lambda x: x.split("_")[0])



#----------------------------------------------
# Pearson correlation of avg ism and avg gradxinp
#----------------------------------------------
corr = df.set_index(['protein', 'dataset']).corr()
plt.figure(figsize=(8,8))
sns.clustermap(corr, annot=True, cmap='coolwarm', fmt=".2f")
plt.savefig("Plots/heatmap_motif_importance_correlation.pdf")
plt.close()


#-----------------------------------
# Pearson correlation between files
#-----------------------------------
def ism_corr_for_4_files(df,ism):
    df_subset_assay=df[["protein","dataset",ism]]
    df_wide=pd.pivot_table(df_subset_assay,index="protein",columns="dataset",values=ism)
    corr=df_wide.corr()
    plt.figure(figsize=(4,4))
    plt.title(ism)
    sns.clustermap(corr, annot=True, cmap='coolwarm', fmt=".2f")
    plt.savefig(f"heatmap_file_correlation_{ism}.pdf")
    plt.close()

ism_corr_for_4_files(df,"ism_cage")
ism_corr_for_4_files(df,"ism_dhs")
ism_corr_for_4_files(df,"ism_starr")
ism_corr_for_4_files(df,"ism_sure")




#----------------------------------------------------------
# Comparison between cell types (color by RE)
#----------------------------------------------------------
def compare_cell_type(df,ism,threshold):
    df_sub=df[["protein","RE","cell_type",ism]]
    df_sub_promoter=df_sub[df_sub["RE"]=="promoters"]
    df_sub_enhancer=df_sub[df_sub["RE"]=="enhancers"]
    # convert to wide format
    df_sub_promoter=df_sub_promoter.pivot(index='protein', columns='cell_type', values=ism).reset_index()
    df_sub_promoter["RE"]="promoters"
    df_sub_enhancer=df_sub_enhancer.pivot(index='protein', columns='cell_type', values=ism).reset_index()
    df_sub_enhancer["RE"]="enhancers"
    df_plot=pd.concat([df_sub_promoter,df_sub_enhancer]).reset_index(drop=True)
    # sns scatter plot
    plt.figure(figsize=(6,6))
    sns.scatterplot(x="hepg2", y="k562", hue="RE", data=df_plot,s=10)
    min_val = min(df_plot[["hepg2", "k562"]].min())
    max_val = max(df_plot[["hepg2", "k562"]].max())
    plt.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--")
    # calculate correlation within promoter, enhancer and for all
    corr_promoter = df_plot[df_plot["RE"] == "promoters"].corr().iloc[0, 1]
    corr_enhancer = df_plot[df_plot["RE"] == "enhancers"].corr().iloc[0, 1]
    corr_all = df_plot.corr().iloc[0, 1]
    plt.text(min_val+(max_val-min_val)*0.3, 
             max_val-(max_val-min_val)*0.1, 
             f"Pearson Correlation:\nPromoter: {corr_promoter:.2f}\nEnhancer: {corr_enhancer:.2f}\nAll: {corr_all:.2f}", fontsize=8)
    # add text annotation for each point
    texts = []
    for i, txt in enumerate(df_plot["protein"]):
        if abs(df_plot["hepg2"][i] - df_plot["k562"][i]) > threshold:
            texts.append(plt.text(df_plot["hepg2"][i], df_plot["k562"][i], txt, fontsize=7))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))
    plt.title(f"HepG2 v.s. K562, {ism}")
    plt.savefig(f"Plots/scatter_hepg2_vs_k562_{ism}.pdf")
    plt.close()

threshold_dict={"ism_cage":0.08,"ism_dhs":0.07,"ism_starr":0.07,"ism_sure":0.25}

for ism in ["ism_cage","ism_dhs","ism_starr","ism_sure"]:
    compare_cell_type(df,ism,threshold_dict[ism])



#----------------------------------------------
# contrast between RE using same track 
# cell line as reproducibility check
#----------------------------------------------

def compare_RE(df,cell_type,ism,threshold1, threshold2):
    df_subset=df[df["dataset"].str.contains(cell_type)]
    df_subset=df_subset[["protein","dataset",ism]]
    # pivot to wide format
    df_wide=pd.pivot_table(df_subset,index="protein",columns="dataset",values=ism)
    plt.figure(figsize=(6,6))
    sns.scatterplot(x=f"promoters_{cell_type}",y=f"enhancers_{cell_type}",data=df_wide,s=10)
    # add diagonal line
    min_val=min(df_wide.min())
    max_val=max(df_wide.max())
    texts = []
    for i, txt in enumerate(df_wide.index):
        diff=df_wide[f"promoters_{cell_type}"][i] - df_wide[f"enhancers_{cell_type}"][i]
        if diff > threshold1 or diff < threshold2:
            texts.append(plt.text(df_wide[f"promoters_{cell_type}"][i], df_wide[f"enhancers_{cell_type}"][i], txt, fontsize=7))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))
    plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
    plt.title(f"promoters vs enhancers {cell_type}: {ism}")
    plt.savefig(f"scatter_promoters_vs_enhancers_{cell_type}_{ism}.pdf")
    plt.close()

threshold1_dict={"ism_cage":0.2,"ism_dhs":0.05,"ism_starr":0.03,"ism_sure":0.3}
threshold2_dict={"ism_cage":-0.05,"ism_dhs":-0.03,"ism_starr":-0.10,"ism_sure":-0.05}

for cell_type in ["hepg2","k562"]:
    for ism in ["ism_cage","ism_dhs","ism_starr","ism_sure"]:
        compare_RE(df,cell_type,ism,threshold1_dict[ism],threshold2_dict[ism])
















# -----------
# plot PCA
# -----------
def plot_pca(df,title,out_name):
    features = ["cage", "dhs", "starr", "sure"]
    x = df.loc[:,features].values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    pca_components = pca.components_
    pc1_label = " ".join(f"{coef:.2f}{feat}" for coef, feat in zip(pca_components[0], features))
    pc2_label = " ".join(f"{coef:.2f}{feat}" for coef, feat in zip(pca_components[1], features))
    
    explained_variance = pca.explained_variance_ratio_
    pc1_explained = f"{explained_variance[0]:.2%}"
    pc2_explained = f"{explained_variance[1]:.2%}"
    
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', data=principalDf)
    for i, txt in enumerate(df['protein']):
        plt.annotate(txt, (principalDf['PC1'][i], principalDf['PC2'][i]), fontsize=8)
    plt.xlabel(f"PC1 = {pc1_label} (EV: {pc1_explained})")
    plt.ylabel(f"PC2 = {pc2_label} (EV: {pc2_explained})")
    # TODO: change suffix
    plt.title(f"{title} classification")
    plt.savefig(f"{out_name}_classification.png")
    plt.close()



plot_pca(df_promoters_hepg2,"Promoters HepG2","promoters_hepg2")
plot_pca(df_promoters_k562,"Promoters K562","promoters_k562")
plot_pca(df_enhancers_hepg2,"Enhancers HepG2","enhancers_hepg2")
plot_pca(df_enhancers_k562,"Enhancers K562","enhancers_k562")

# nohup python3 quantify_motif_effect.py > quantify_motif_effect.out &