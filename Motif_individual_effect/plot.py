import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from adjustText import adjust_text
from loguru import logger

#----------------------------------------------
# Helper functions
#----------------------------------------------



#----------------------------------------------
# Pearson correlation of avg ism and avg gradxinp, dstat ism and dstat gradxinp
#----------------------------------------------
def generate_heatmap(df, column_filters, output_file):
    columns = df.columns[df.columns.str.contains(column_filters[0]) | df.columns.str.contains(column_filters[1])]
    df_sub = df[["protein", "dataset"] + list(columns)]
    corr = df_sub.set_index(['protein', 'dataset']).corr()
    plt.figure(figsize=(8, 8))
    sns.clustermap(corr, annot=True, cmap='coolwarm', fmt=".2f")
    plt.savefig(output_file)
    plt.close()


# Method 1: agregate by file
df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_file.csv")
generate_heatmap(df, ["avg_ism", "avg_gradxinp"], "Plots/heatmap_motif_importance_correlation_aggregate_by_file_avg_ism_vs_avg_gradxinp.pdf")
generate_heatmap(df, ["dstat_ism", "dstat_gradxinp"], "Plots/heatmap_motif_importance_correlation_aggregate_by_file_dstat_ism_vs_dstat_gradxinp.pdf")
generate_heatmap(df, ["dstat_ism", "avg_ism"], "Plots/heatmap_motif_importance_correlation_aggregate_by_file_dstat_ism_vs_avg_ism.pdf")
generate_heatmap(df, ["dstat_gradxinp", "avg_gradxinp"], "Plots/heatmap_motif_importance_correlation_aggregate_by_file_dstat_gradxinp_vs_avg_gradxinp.pdf")

# Method 2: agregate by cell type
df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_cell_type.csv")
generate_heatmap(df, ["avg_ism", "avg_gradxinp"], "Plots/heatmap_motif_importance_correlation_aggregate_by_cell_avg_ism_vs_avg_gradxinp.pdf")
generate_heatmap(df, ["dstat_ism", "dstat_gradxinp"], "Plots/heatmap_motif_importance_correlation_aggregate_by_cell_dstat_ism_vs_dstat_gradxinp.pdf")
generate_heatmap(df, ["dstat_ism", "avg_ism"], "Plots/heatmap_motif_importance_correlation_aggregate_by_cell_dstat_ism_vs_avg_ism.pdf")
generate_heatmap(df, ["dstat_gradxinp", "avg_gradxinp"], "Plots/heatmap_motif_importance_correlation_aggregate_by_cell_dstat_gradxinp_vs_avg_gradxinp.pdf")



#-----------------------------------
# Pearson correlation between files
#-----------------------------------


def ism_corr_for_multiple_files(df, isms, out_suffix):
    _, axes = plt.subplots(2, 2, figsize=(12,12)) 
    axes = axes.flatten()  # Flatten the array for easier indexing
    for i, ism in enumerate(isms):
        df_subset_assay = df[["protein", "dataset", ism]]
        df_wide = pd.pivot_table(df_subset_assay, index="protein", columns="dataset", values=ism)
        corr = df_wide.corr()
        sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", ax=axes[i])
        axes[i].set_title(ism)
    plt.tight_layout()
    plt.savefig(f"heatmap_file_correlations_{out_suffix}.pdf")
    plt.close()

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_file.csv")
df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["RE"]=df["dataset"].apply(lambda x: x.split("_")[0])

ism_corr_for_multiple_files(df, ["avg_ism_cage", "avg_gradxinp_cage", "dstat_ism_cage", "dstat_gradxinp_cage"],"cage")
ism_corr_for_multiple_files(df, ["avg_ism_dhs", "avg_gradxinp_dhs", "dstat_ism_dhs", "dstat_gradxinp_dhs"],"dhs")
ism_corr_for_multiple_files(df, ["avg_ism_starr", "avg_gradxinp_starr", "dstat_ism_starr", "dstat_gradxinp_starr"],"starr")
ism_corr_for_multiple_files(df, ["avg_ism_sure", "avg_gradxinp_sure", "dstat_ism_sure", "dstat_gradxinp_sure"],"sure")



#----------------------------------------------------------
# Comparison between cell types (color by RE)
#----------------------------------------------------------
def compare_between_cell_type_separate_by_re(df,ism,threshold,mode):
    df_sub = df[["protein", "RE", "cell_type", ism]]
    df_sub_promoter = df_sub[df_sub["RE"] == "promoters"]
    df_sub_enhancer = df_sub[df_sub["RE"] == "enhancers"]
    df_sub_promoter = df_sub_promoter.pivot(index='protein', columns='cell_type', values=ism).reset_index()
    df_sub_promoter["RE"] = "promoters"
    df_sub_enhancer = df_sub_enhancer.pivot(index='protein', columns='cell_type', values=ism).reset_index()
    df_sub_enhancer["RE"] = "enhancers"
    df_plot = pd.concat([df_sub_promoter, df_sub_enhancer]).reset_index(drop=True)
    # replace nan with 0
    df_plot.fillna(0, inplace=True)
    plt.figure(figsize=(6, 6))
    max_difference = max(abs(df_plot["hepg2"] - df_plot["k562"]))
    # Plot each dot with individual color and transparency
    for i in range(len(df_plot)):
        difference = abs(df_plot["hepg2"][i] - df_plot["k562"][i])
        alpha_value = difference / max_difference # Calculate transparency
        color = '#1f77b4' if df_plot["RE"][i] == "promoters" else '#ff7f0e'
        plt.scatter(df_plot["hepg2"][i], df_plot["k562"][i], color=color, alpha=alpha_value, s=3)
    min_val = min(df_plot[["hepg2", "k562"]].min())
    max_val = max(df_plot[["hepg2", "k562"]].max())
    plt.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--")
    # Calculate correlations
    corr_promoter = df_plot[df_plot["RE"] == "promoters"].corr().iloc[0, 1]
    corr_enhancer = df_plot[df_plot["RE"] == "enhancers"].corr().iloc[0, 1]
    corr_all = df_plot.corr().iloc[0, 1]
    # Add correlation text
    plt.text(min_val + (max_val - min_val) * 0.3, 
             max_val - (max_val - min_val) * 0.1, 
             f"Pearson Correlation:\nPromoter: {corr_promoter:.2f}\nEnhancer: {corr_enhancer:.2f}\nAll: {corr_all:.2f}", fontsize=8)
    # Add text annotations
    texts = []
    for i, txt in enumerate(df_plot["protein"]):
        if mode == "avg_ism":
            if abs(df_plot["hepg2"][i] - df_plot["k562"][i]) > threshold:
                texts.append(plt.text(df_plot["hepg2"][i], df_plot["k562"][i], txt, fontsize=4))
        elif mode == "dstat_ism":
            condition1 = (df_plot["hepg2"][i] * df_plot["k562"][i] < 0) and (abs(df_plot["hepg2"][i] - df_plot["k562"][i]) > 0.4)
            condition2 = (df_plot["hepg2"][i] + df_plot["k562"][i] > 0.9)
            condition3 = (df_plot["hepg2"][i] + df_plot["k562"][i] < -0.9)
            if condition1 or condition2 or condition3:
                texts.append(plt.text(df_plot["hepg2"][i], df_plot["k562"][i], txt, fontsize=4))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5))
    # Plot title and save
    plt.title(f"HepG2 vs. K562, {ism}")
    plt.savefig(f"Plots/scatter_hepg2_vs_k562_{ism}.pdf")
    plt.close()



df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_file.csv")
df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["RE"]=df["dataset"].apply(lambda x: x.split("_")[0])
threshold_dict={"avg_ism_cage":0.09,"avg_ism_dhs":0.07,"avg_ism_starr":0.07,"avg_ism_sure":0.25}
for ism in threshold_dict.keys():
    compare_between_cell_type_separate_by_re(df,ism,threshold_dict[ism],mode="avg_ism")
for ism in ["dstat_ism_cage","dstat_ism_dhs","dstat_ism_starr","dstat_ism_sure"]:
    compare_between_cell_type_separate_by_re(df,ism,None,mode="dstat_ism")





#----------------------------------------------------------
# Comparison between cell types (Don't separate by RE)
#----------------------------------------------------------
def compare_between_cell_types(df,ism,threshold,mode):
    df_sub = df[["protein","dataset", ism]]
    df_plot = df_sub.pivot(index='protein', columns='dataset', values=ism).reset_index()
    # replace nan with 0
    df_plot.fillna(0, inplace=True)
    plt.figure(figsize=(6, 6))
    max_difference = max(abs(df_plot["hepg2"] - df_plot["k562"]))
    # Plot each dot with individual color and transparency
    for i in range(len(df_plot)):
        difference = df_plot["hepg2"][i] - df_plot["k562"][i]
        if difference > 0:
            color = '#1f77b4' 
        elif difference < 0:
            color = '#ff7f0e'
        difference = abs(difference)
        alpha_value = difference / max_difference 
        plt.scatter(df_plot["hepg2"][i], df_plot["k562"][i], color=color, alpha=alpha_value, s=5)
    # Add text
    texts = []
    for i, txt in enumerate(df_plot["protein"]):
        if mode == "avg_ism":
            if abs(df_plot["hepg2"][i] - df_plot["k562"][i]) > threshold:
                texts.append(plt.text(df_plot["hepg2"][i], df_plot["k562"][i], txt, fontsize=6))
        elif mode == "dstat_ism":
            condition1 = (df_plot["hepg2"][i] * df_plot["k562"][i] < 0) and (abs(df_plot["hepg2"][i] - df_plot["k562"][i]) > 0.4)
            condition2 = (df_plot["hepg2"][i] + df_plot["k562"][i] > 0.9)
            condition3 = (df_plot["hepg2"][i] + df_plot["k562"][i] < -0.9)
            if condition1 or condition2 or condition3:
                texts.append(plt.text(df_plot["hepg2"][i], df_plot["k562"][i], txt, fontsize=6))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5))
    # Plot title and save
    min_val = min(df_plot[["hepg2", "k562"]].min())
    max_val = max(df_plot[["hepg2", "k562"]].max())
    plt.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--")
    # add vertical line
    plt.axvline(x=0, color='black')
    # add horizontal line
    plt.axhline(y=0, color='black')
    plt.title(f"HepG2 vs. K562, {ism}")
    plt.xlabel("Motif_importance_for_HepG2")
    plt.ylabel("Motif_importance_for_K562")
    plt.savefig(f"Plots/scatter_hepg2_vs_k562_no_re_split_{ism}.pdf")
    plt.close()


df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_cell_type.csv")
threshold_dict={"avg_ism_cage":0.06,"avg_ism_dhs":0.04,"avg_ism_starr":0.04,"avg_ism_sure":0.20}
for ism in threshold_dict.keys():
    compare_between_cell_types(df,ism,threshold_dict[ism],mode="avg_ism")
for ism in ["dstat_ism_cage","dstat_ism_dhs","dstat_ism_starr","dstat_ism_sure"]:
    compare_between_cell_types(df,ism,None,mode="dstat_ism")






#----------------------------------------------
# contrast ism between RE using same track 
# cell line as reproducibility check
#----------------------------------------------

def compare_RE_ism(df,cell_type,ism,threshold1, threshold2):
    df_subset=df[df["dataset"].str.contains(cell_type)]
    df_subset=df_subset[["protein","dataset",ism]]
    # pivot to wide format
    df_wide=pd.pivot_table(df_subset,index="protein",columns="dataset",values=ism)
    plt.figure(figsize=(6,6))
    sns.scatterplot(x=f"promoters_{cell_type}",y=f"enhancers_{cell_type}",data=df_wide,s=10)
    texts = []
    for i, txt in enumerate(df_wide.index):
        diff=df_wide[f"promoters_{cell_type}"][i] - df_wide[f"enhancers_{cell_type}"][i]
        if diff > threshold1 or diff < threshold2:
            texts.append(plt.text(df_wide[f"promoters_{cell_type}"][i], df_wide[f"enhancers_{cell_type}"][i], txt, fontsize=5))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))
    # add diagonal line
    min_val=min(df_wide.min())
    max_val=max(df_wide.max())
    plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
    plt.title(f"promoters vs enhancers {cell_type}: {ism}")
    plt.savefig(f"Plots/scatter_promoters_vs_enhancers_{cell_type}_{ism}.pdf")
    plt.close()

threshold1_dict={"avg_ism_cage":0.3,"avg_ism_dhs":0.05,"avg_ism_starr":0.03,"avg_ism_sure":0.2}
threshold2_dict={"avg_ism_cage":-0.05,"avg_ism_dhs":-0.03,"avg_ism_starr":-0.10,"avg_ism_sure":-0.2}

for cell_type in ["hepg2","k562"]:
    for ism in ["avg_ism_sure"]:
        compare_RE_ism(df,cell_type,ism,threshold1_dict[ism],threshold2_dict[ism])





#----------------------------------------------
# contrast dstat between RE using same track 
# cell line as reproducibility check
#----------------------------------------------

def scatter_plot_alpha_by_diff(df_plot,x_colname,y_colname,ism,annot_condition_colname="annotate"):
    max_difference = max(abs(df_plot[x_colname] - df_plot[y_colname]))
    plt.figure(figsize=(6, 6))
    # Plot each dot with individual color and transparency
    for i in range(df_plot.shape[0]):
        difference = abs(df_plot[x_colname][i] - df_plot[y_colname][i])
        alpha_value = difference / max_difference
        plt.scatter(df_plot[x_colname][i], df_plot[y_colname][i], color='#1f77b4', alpha=alpha_value, s=3)
    # Add text annotations
    texts = []
    df_annotate = df_plot[df_plot[annot_condition_colname]==True].reset_index(drop=True)
    for i, txt in enumerate(df_annotate["protein"]):
        texts.append(plt.text(df_annotate[x_colname][i], df_annotate[y_colname][i], txt, fontsize=4))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5))
    # add diagonal line
    min_val = min(df_plot[[x_colname, y_colname]].min())
    max_val = max(df_plot[[x_colname, y_colname]].max())
    # add vertical line
    plt.axvline(x=0, color='black')
    # add horizontal line
    plt.axhline(y=0, color='black')
    plt.xlabel(f"{x_colname}")
    plt.ylabel(f"{y_colname}")
    plt.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--")
    plt.title(f"{x_colname} vs {y_colname}, {ism}")
    plt.savefig(f"Plots/scatter_{x_colname}_vs_{y_colname}_{ism}.pdf")
    plt.close()



def compare_RE_dstat(df,cell_type,ism):
    x_colname=f"promoters_{cell_type}"
    y_colname=f"enhancers_{cell_type}"
    df_subset=df[df["dataset"].str.contains(cell_type)]
    df_subset=df_subset[["protein","dataset",ism]]
    # pivot to wide format
    df_wide=pd.pivot_table(df_subset,index="protein",columns="dataset",values=ism)
    df_wide["protein"]=df_wide.index
    df_wide.fillna(0,inplace=True)
    # calculate whether to annotate the row
    condition1 = (df_wide[x_colname] * df_wide[y_colname] < 0) & (abs(df_wide[x_colname] - df_wide[y_colname]) > 0.4)
    condition2 = (df_wide[x_colname] + df_wide[y_colname] > 0.9)
    condition3 = (df_wide[x_colname] + df_wide[y_colname] < -0.9)
    df_wide["annotate"]=condition1 | condition2 | condition3
    scatter_plot_alpha_by_diff(df_wide,x_colname,y_colname,ism)
    


df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_file.csv")
df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["RE"]=df["dataset"].apply(lambda x: x.split("_")[0])
for cell_type in ["hepg2","k562"]:
    for ism in ["dstat_ism_cage","dstat_ism_dhs","dstat_ism_starr","dstat_ism_sure"]:
        compare_RE_dstat(df,cell_type,ism)







#----------------------------------------------
# contrast between tracks 
# cell line as reproducibility check
#----------------------------------------------


df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_cell_type.csv")
# select only HepG2
for cell_type in ["hepg2","k562"]:
    df_sub=df[df["dataset"].str.contains(cell_type)].reset_index(drop=True)
    assays=["cage","dhs","starr","sure"]
    for i in range(3):
        for j in range(i+1,4):
            x_colname=f"dstat_ism_{assays[i]}"
            y_colname=f"dstat_ism_{assays[j]}"
        condition1 = (df_sub[x_colname] * df_sub[y_colname] < 0) & (abs(df_sub[x_colname] - df_sub[y_colname]) > 0.4)
        condition2 = (df_sub[x_colname] + df_sub[y_colname] > 0.9)
        condition3 = (df_sub[x_colname] + df_sub[y_colname] < -0.9)
        df_sub["annotate"]=condition1 | condition2 | condition3
        scatter_plot_alpha_by_diff(df_sub,x_colname,y_colname,cell_type)










# -----------
# plot PCA
# -----------
def plot_pca(df,file,method):
    df=df[df["dataset"]==file]
    features=["cage","dhs","starr","sure"]
    features_prefixed=[f"{method}_{x}" for x in features]
    x = df.loc[:,features_prefixed].values
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    pca_components = pca.components_
    pc1_label = " ".join(f"{coef:.2f}{feat}" for coef, feat in zip(pca_components[0], features))
    pc2_label = " ".join(f"{coef:.2f}{feat}" for coef, feat in zip(pca_components[1], features))
    explained_variance = pca.explained_variance_ratio_
    pc1_explained = f"{explained_variance[0]:.2%}"
    pc2_explained = f"{explained_variance[1]:.2%}"
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    plt.figure(figsize=(6, 6))
    sns.scatterplot(x='PC1', y='PC2', data=principalDf,s=3)
    for i, txt in enumerate(df['protein']):
        plt.annotate(txt, (principalDf['PC1'][i], principalDf['PC2'][i]), fontsize=4)
    plt.xlabel(f"PC1 = {pc1_label} (EV: {pc1_explained})")
    plt.ylabel(f"PC2 = {pc2_label} (EV: {pc2_explained})")
    plt.title(f"PCA of {file} {method}")
    plt.savefig(f"pca_{file}_{method}.pdf")
    plt.close()

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_file.csv")
plot_pca(df,"promoters_hepg2","avg_ism")
plot_pca(df,"promoters_hepg2","dstat_ism")
plot_pca(df,"enhancers_hepg2","avg_ism")
plot_pca(df,"enhancers_hepg2","dstat_ism")
plot_pca(df,"promoters_k562","avg_ism")
plot_pca(df,"promoters_k562","dstat_ism")
plot_pca(df,"enhancers_k562","avg_ism")
plot_pca(df,"enhancers_k562","dstat_ism")


df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_cell_type.csv")
plot_pca(df,"hepg2","dstat_ism")
plot_pca(df,"k562","dstat_ism")


