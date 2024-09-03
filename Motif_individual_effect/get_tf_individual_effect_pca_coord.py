from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt





def read_tf_individual_effect(by):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_individual_effect_by_{by}.csv")
    columns=[col for col in df.columns if "dstat_ism" in col]
    df=df[columns+["protein","dataset"]]
    return df



def write_pca_coord_by_cell(cell_type):
    df=read_tf_individual_effect("cell_type")
    df=df[df["dataset"]==cell_type].reset_index(drop=True)
    df["tf_type"]=np.where(df["protein"].isin(sub_tfs),"sub",np.where(df["protein"].isin(super_tfs),"super","other"))
    pca = PCA(n_components=2)
    features=["dstat_ism_"+x for x in ["cage","dhs","starr","sure"]]
    x = df.loc[:,features].values
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    principalDf["tf_type"]=df["tf_type"]
    principalDf["protein"]=df["protein"]
    principalDf.to_csv(f"tf_pca_coord_with_sub_super_{cell_type}.csv",index=False)
    


def write_pca_coord_by_file(file_name):
    df=read_tf_individual_effect("file")
    df=df[df["dataset"]==file_name].reset_index(drop=True)
    sub_tfs,super_tfs=get_sub_super_annotation([file_name])
    df["tf_type"]=np.where(df["protein"].isin(sub_tfs),"sub",np.where(df["protein"].isin(super_tfs),"super","other"))
    pca = PCA(n_components=2)
    features=["dstat_ism_"+x for x in ["cage","dhs","starr","sure"]]
    x = df.loc[:,features].values
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    principalDf["tf_type"]=df["tf_type"]
    principalDf["protein"]=df["protein"]
    principalDf.to_csv(f"tf_pca_coord_with_sub_super_{file_name}.csv",index=False)



write_pca_coord_by_cell("hepg2")
write_pca_coord_by_cell("k562")
write_pca_coord_by_file("enhancers_hepg2")
write_pca_coord_by_file("promoters_hepg2")
write_pca_coord_by_file("enhancers_k562")
write_pca_coord_by_file("promoters_k562")



cell="k562"
df=read_tf_individual_effect("file_chip_true")
# remove rows with nan
df=df.dropna()
df["cell"]=df["dataset"].apply(lambda x: x.split("_")[1])
df["re"]=df["dataset"].apply(lambda x: x.split("_")[0])
# subset to contain only hepg2 cell
df=df[df["cell"]==cell].reset_index(drop=True)
# split df to df_promoter and df_enhancer
df_promoter=df[df["re"]=="promoters"].reset_index(drop=True)
df_enhancer=df[df["re"]=="enhancers"].reset_index(drop=True)
# rename columns
df_promoter.columns=[col.replace("dstat_ism","promoter") for col in df_promoter.columns]
df_enhancer.columns=[col.replace("dstat_ism","enhancer") for col in df_enhancer.columns]
# remover columns 'dataset', 'cell', 're', set indext to 'protein'
df_promoter=df_promoter.drop(columns=["dataset","cell","re"]).set_index("protein")
df_enhancer=df_enhancer.drop(columns=["dataset","cell","re"]).set_index("protein")

# merge df_promoter and df_enhancer
df=pd.merge(df_promoter,df_enhancer,left_index=True,right_index=True,how="inner")
# perform PCA
pca = PCA(n_components=2)
x = df.values
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
principalDf["protein"]=df.index
principalDf.to_csv(f"tf_pca_coord_{cell}_chip_true.csv",index=False)
# get meaning of pc1 and pc2
print(pca.components_)
# get explained variance ratio
print(pca.explained_variance_ratio_)


plt.figure(figsize=(6, 6))
sns.scatterplot(x='PC1', y='PC2', data=principalDf,s=3)
for i, txt in enumerate(principalDf['protein']):
    plt.annotate(txt, (principalDf['PC1'][i], principalDf['PC2'][i]), fontsize=4)

plt.savefig(f"pca_{cell}.pdf")
plt.close()


