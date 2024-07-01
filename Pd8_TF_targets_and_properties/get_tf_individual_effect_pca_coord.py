from sklearn.decomposition import PCA
import pandas as pd
import numpy as np




def read_tf_individual_effect(by):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_{by}.csv")
    columns=[col for col in df.columns if "dstat_ism" in col]
    df=df[columns+["protein","dataset"]]
    return df



def get_sub_super_annotation(datasets):
    df_tf_property=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_additivity_property.csv",index_col=0)
    df_tf_property=df_tf_property[datasets]
    df_tf_property["sub_counts"]=df_tf_property.apply(lambda x: x.value_counts().get("sub",0),axis=1)
    df_tf_property["super_counts"]=df_tf_property.apply(lambda x: x.value_counts().get("super",0),axis=1)
    sub_tfs =  df_tf_property[(df_tf_property["sub_counts"]>0) & (df_tf_property["super_counts"]==0)].index.to_list()
    super_tfs = df_tf_property[(df_tf_property["super_counts"]>0) & (df_tf_property["sub_counts"]==0)].index.to_list()
    return sub_tfs,super_tfs


def write_pca_coord_by_cell(cell_type):
    df=read_tf_individual_effect("cell_type")
    df=df[df["dataset"]==cell_type].reset_index(drop=True)
    sub_tfs,super_tfs=get_sub_super_annotation([f"enhancers_{cell_type}", f"promoters_{cell_type}"])
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