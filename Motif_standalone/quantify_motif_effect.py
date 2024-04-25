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
    x = df.loc[:,["cage", "dhs", "starr", "sure"]].values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', data=principalDf)
    for i, txt in enumerate(df['protein']):
        plt.annotate(txt, (principalDf['PC1'][i], principalDf['PC2'][i]), fontsize=8)
    plt.title(title)
    plt.savefig(out_name)
    plt.close()



df_promoters_hepg2=read_data("promoters_hepg2", [0,2,4,6])
df_promoters_k562=read_data("promoters_k562", [1,3,5,7])
df_enhancers_hepg2=read_data("enhancers_hepg2", [0,2,4,6])
df_enhancers_k562=read_data("enhancers_k562", [1,3,5,7])


df_promoters_hepg2=calc_avg_ism(df_promoters_hepg2)
df_promoters_k562=calc_avg_ism(df_promoters_k562)
df_enhancers_hepg2=calc_avg_ism(df_enhancers_hepg2)
df_enhancers_k562=calc_avg_ism(df_enhancers_k562)

plot_pca(df_promoters_hepg2,"Promoters HepG2","promoters_hepg2.png")
plot_pca(df_promoters_k562,"Promoters K562","promoters_k562.png")
plot_pca(df_enhancers_hepg2,"Enhancers HepG2","enhancers_hepg2.png")
plot_pca(df_enhancers_k562,"Enhancers K562","enhancers_k562.png")