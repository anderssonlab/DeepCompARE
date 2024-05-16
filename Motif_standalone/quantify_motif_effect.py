import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


prefix="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_"

def read_data(re_type, track_list):
    df_list=[]
    for track in track_list:
        df=pd.read_csv(f"{prefix}{re_type}_track{track}.csv")
        df["dataset"]=f"track{track}"
        df_list.append(df)
    df_res=pd.concat(df_list)
    df_res=df_res[df_res["chip_evidence"]==True]
    return df_res


def calc_avg_ism(df):
    df_res=df.groupby(["dataset","protein"])["ism_motif"].mean().reset_index()
    df_res=df_res.pivot(index='protein', columns='dataset', values='ism_motif').reset_index()
    df_res.columns=["protein", "cage", "dhs", "starr", "sure"]
    return df_res


# perform A PCA on this df_promoters_ism, plot the PCA
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



df_promoters_hepg2=read_data("promoters_hepg2", [8,10,12,14])
df_promoters_k562=read_data("promoters_k562", [9,11,13,15])
df_enhancers_hepg2=read_data("enhancers_hepg2", [8,10,12,14])
df_enhancers_k562=read_data("enhancers_k562", [9,11,13,15])


df_promoters_hepg2=calc_avg_ism(df_promoters_hepg2)
df_promoters_k562=calc_avg_ism(df_promoters_k562)
df_enhancers_hepg2=calc_avg_ism(df_enhancers_hepg2)
df_enhancers_k562=calc_avg_ism(df_enhancers_k562)

plot_pca(df_promoters_hepg2,"Promoters HepG2","promoters_hepg2")
plot_pca(df_promoters_k562,"Promoters K562","promoters_k562")
plot_pca(df_enhancers_hepg2,"Enhancers HepG2","enhancers_hepg2")
plot_pca(df_enhancers_k562,"Enhancers K562","enhancers_k562")

# nohup python3 quantify_motif_effect.py > quantify_motif_effect.out &