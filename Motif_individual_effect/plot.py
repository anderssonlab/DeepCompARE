import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from adjustText import adjust_text
from loguru import logger



def scatter_plot_alpha_by_diff(df_plot,x_colname,y_colname,ism,annot_condition_colname="annotate"):
    max_difference = max(abs(df_plot[x_colname] - df_plot[y_colname]))
    plt.figure(figsize=(6, 6))
    # Plot each dot with individual color and transparency
    for i in range(df_plot.shape[0]):
        if df_plot[x_colname][i] - df_plot[y_colname][i]>0:
            color='#1f77b4'
        else:
            color='#ff7f0e'
        difference = abs(df_plot[x_colname][i] - df_plot[y_colname][i])
        alpha_value = min((1.5* difference / max_difference),1)
        plt.scatter(df_plot[x_colname][i], df_plot[y_colname][i], color=color, alpha=alpha_value, s=5)
    # Add text annotations
    texts = []
    df_annotate = df_plot[df_plot[annot_condition_colname]==True].reset_index(drop=True)
    for i, txt in enumerate(df_annotate["protein"]):
        texts.append(plt.text(df_annotate[x_colname][i], df_annotate[y_colname][i], txt, fontsize=5))
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




#----------------------------------------------
# Pearson correlation of avg ism and dstat ism
#----------------------------------------------
def generate_heatmap(df,output_file):
    # TODO: choose dstat_ism or avg_ism
    # avg can cluster by in vivo and in vitro, by cell v.s. by file performs similarly
    columns = df.columns[df.columns.str.contains("avg_ism")]
    df_sub = df[["protein", "dataset"] + list(columns)]
    corr = df_sub.set_index(['protein', 'dataset']).corr()
    sns.clustermap(corr, annot=True, cmap='coolwarm', fmt=".2f",figsize=(6, 6))
    plt.savefig(output_file)
    plt.close()


# Method 1: agregate by cell type
df=pd.read_csv("tf_individual_effect_by_cell_type.csv")
generate_heatmap(df, "Plots/heatmap_motif_importance_correlation_avg_by_cell.pdf")

# Method 2: agregate by file
df=pd.read_csv("tf_individual_effect_by_file.csv")
generate_heatmap(df,"Plots/heatmap_motif_importance_correlation_avg_by_file.pdf")

# final decision: avg, by cell

#-----------------------------------
# Pearson correlation between files
#-----------------------------------


def ism_corr_for_multiple_files(df, ism, out_suffix):
    df_subset_assay = df[["protein", "dataset", ism]]
    df_wide = pd.pivot_table(df_subset_assay, index="protein", columns="dataset", values=ism)
    corr = df_wide.corr()
    clustermap = sns.clustermap(corr, annot=True, cmap='coolwarm', fmt=".2f", figsize=(6, 6))
    clustermap.fig.suptitle(f"Correlation of {ism} effect", y=1.05)
    clustermap.savefig(f"Plots/heatmap_file_correlations_{out_suffix}.pdf")
    plt.close()

df=pd.read_csv("tf_individual_effect_by_file.csv")
df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["RE"]=df["dataset"].apply(lambda x: x.split("_")[0])

ism_corr_for_multiple_files(df, "avg_ism_cage","cage")
ism_corr_for_multiple_files(df, "avg_ism_dhs","dhs")
ism_corr_for_multiple_files(df, "avg_ism_starr","starr")
ism_corr_for_multiple_files(df, "avg_ism_sure","sure")

#----------------------------------------------------------
# Comparison between cell types (Don't separate by RE)
#----------------------------------------------------------
def compare_between_cell_types(df,ism):
    df_sub = df[["protein","dataset", ism]]
    df_plot = df_sub.pivot(index='protein', columns='dataset', values=ism).reset_index()
    condition1 = (df_plot["hepg2"] * df_plot["k562"] < 0) & (abs(df_plot["hepg2"] - df_plot["k562"]) > 0.3)
    condition2 = (df_plot["hepg2"]>0.5) & (df_plot["k562"] > 0.5)
    condition3 = (df_plot["hepg2"] < -0.5) & (df_plot["k562"] < -0.5)
    df_plot["annotate"]= condition1 | condition2 | condition3
    scatter_plot_alpha_by_diff(df_plot,"hepg2","k562",ism)


df=pd.read_csv("tf_individual_effect_by_cell_type.csv")
for ism in ["dstat_ism_cage","dstat_ism_dhs","dstat_ism_starr","dstat_ism_sure"]:
    compare_between_cell_types(df,ism)




#----------------------------------------------
# contrast between tracks 
# cell line as reproducibility check
#----------------------------------------------

df=pd.read_csv("tf_individual_effect_by_cell_type.csv")
for cell_type in ["hepg2","k562"]:
    df_sub=df[df["dataset"].str.contains(cell_type)].reset_index(drop=True)
    assays=["cage","dhs","starr","sure"]
    for i in range(3):
        for j in range(i+1,4):
            x_colname=f"dstat_ism_{assays[i]}"
            y_colname=f"dstat_ism_{assays[j]}"
            logger.info(f"Plotting {x_colname} vs {y_colname} for {cell_type}")
            condition1 = (abs(df_sub[x_colname] - df_sub[y_colname]) > 0.15)
            condition2 = (df_sub[x_colname]>0.5) & (df_sub[y_colname] > 0.5)
            condition3 = (df_sub[x_colname]<-0.5) & (df_sub[y_colname] < -0.5)
            df_sub["annotate"]=condition1 | condition2 | condition3
            scatter_plot_alpha_by_diff(df_sub,x_colname,y_colname,cell_type)









#----------------------------------------------
# contrast dstat between RE using same track 
# cell line as reproducibility check
#----------------------------------------------


def compare_RE(df,cell_type,ism):
    x_colname=f"promoters_{cell_type}"
    y_colname=f"enhancers_{cell_type}"
    df_subset=df[df["dataset"].str.contains(cell_type)]
    df_subset=df_subset[["protein","dataset",ism]]
    # pivot to wide format
    df_wide=pd.pivot_table(df_subset,index="protein",columns="dataset",values=ism)
    df_wide["protein"]=df_wide.index
    df_wide.fillna(0,inplace=True)
    # calculate whether to annotate the row
    condition1 = (abs(df_wide[x_colname] - df_wide[y_colname]) > 0.2)
    condition2 = (df_wide[x_colname]>0.5) & (df_wide[y_colname] > 0.5)
    condition3 = (df_wide[x_colname]<-0.5) & (df_wide[y_colname] < -0.5)
    df_wide["annotate"]= condition1 | condition2 | condition3
    scatter_plot_alpha_by_diff(df_wide,x_colname,y_colname,ism)
    


df=pd.read_csv("tf_individual_effect_by_file.csv")
df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["RE"]=df["dataset"].apply(lambda x: x.split("_")[0])
for cell_type in ["hepg2","k562"]:
    for ism in ["dstat_ism_cage","dstat_ism_dhs","dstat_ism_starr","dstat_ism_sure"]:
        compare_RE(df,cell_type,ism)








# -----------
# plot PCA
# -----------
def plot_pca(df):
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






df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_pca_coord_hepg2_chip_true.csv")
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





