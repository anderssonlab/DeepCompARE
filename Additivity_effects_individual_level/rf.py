import pandas as pd
import numpy as np
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import split_dimer
from stat_tests import bin_and_label


#----------------------------
# Helper functions
#----------------------------

def get_sub_super_tfs(suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_pca_coord_with_sub_super_{suffix}.csv")
    sub_tfs=df[df["tf_type"]=="sub"]["protein"].to_list()
    super_tfs=df[df["tf_type"]=="super"]["protein"].to_list()
    sub_tfs=split_dimer(sub_tfs)
    super_tfs=split_dimer(super_tfs)
    intersection=set(sub_tfs).intersection(set(super_tfs))
    sub_tfs=[tf for tf in sub_tfs if tf not in intersection]
    super_tfs=[tf for tf in super_tfs if tf not in intersection]
    return sub_tfs,super_tfs

def add_tf_annotation(df,sub_tfs,super_tfs):
    # Add TF annotation to df
    df['sub-additive']=df["motif"].isin(sub_tfs).astype(int)
    df['super-additive']=df["motif"].isin(super_tfs).astype(int)
    df=df.join(df_dispersion.set_index("gene"),on="motif",how="left")
    df.dropna(subset=["gini"],inplace=True)
    return df


def read_tfbs_maf(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_{file_suffix}.csv",header=None)
    df.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]
    df["dataset"]=file_suffix
    df["log10_AF"]=np.log10(df["AF"])
    return df



def read_maf_effect_size(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Effect_size_vs_maf/maf_with_effect_size_{file_suffix}.csv",header=None,index_col=0)
    df.reset_index(drop=True,inplace=True)
    df.columns=["Chromosome","Start","End","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
    return df



def join_tfbs_maf_effect_size(file_suffix):
    tfbs_maf=read_tfbs_maf(file_suffix)
    maf_effect_size=read_maf_effect_size(file_suffix)
    df=tfbs_maf.join(maf_effect_size.set_index(["Chromosome","Start","End","ID","REF","ALT","AF"]),on=["Chromosome","Start","End","ID","REF","ALT","AF"],how="inner")
    return df



#----------------------------
# Read in data
#----------------------------

# dispersion data
df_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.hepG2.tab",sep="\t")
df_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/TFs.dispersionEstimates.k562.tab",sep="\t")
df_dispersion=pd.concat([df_dispersion_hepg2,df_dispersion_k562],axis=0).reset_index(drop=True)
df_dispersion=df_dispersion.drop_duplicates().reset_index(drop=True)
df_dispersion=df_dispersion[["gene","gini"]]

# df_tfbs_maf_es is for effect size analysis
df_tfbs_maf_es=pd.concat([join_tfbs_maf_effect_size("enhancers_hepg2"),
                          join_tfbs_maf_effect_size("promoters_hepg2"),
                          join_tfbs_maf_effect_size("enhancers_k562"),
                          join_tfbs_maf_effect_size("promoters_k562")],
                          axis=0).reset_index(drop=True)
df_tfbs_maf_es=df_tfbs_maf_es[df_tfbs_maf_es["chip_evidence"]==True].reset_index(drop=True)
df_tfbs_maf_es=bin_and_label(df_tfbs_maf_es,"log10_AF",[-np.inf,-4,-2,0])

#----------------------------
# Add annotation
#----------------------------
df_tfbs_maf_es["promoter"]=df_tfbs_maf_es["dataset"].str.contains("promoters").astype(int)
df_tfbs_maf_es["enhancer"]=df_tfbs_maf_es["dataset"].str.contains("enhancers").astype(int)
df_tfbs_maf_es["hepg2"]=df_tfbs_maf_es["dataset"].str.contains("hepg2").astype(int)
df_tfbs_maf_es["k562"]=df_tfbs_maf_es["dataset"].str.contains("k562").astype(int)


# deduplicate (focus on integrating functional annotations)
dup_mask=df_tfbs_maf_es.duplicated(subset=["Chromosome","Start","End","ID","REF","ALT"],keep=False)
df_tfbs_maf_es_nodup=df_tfbs_maf_es[~dup_mask].reset_index(drop=True)
df_tfbs_maf_es_dup=df_tfbs_maf_es[dup_mask].reset_index(drop=True)
df_tfbs_maf_es_dup_annotation=df_tfbs_maf_es_dup.groupby(["Chromosome",
                                                          "Start",
                                                          "End",
                                                          "ID",
                                                          "REF",
                                                          "ALT",
                                                          "AF"]).agg({"promoter":"max","enhancer":"max","hepg2":"max","k562":"max"}).reset_index()
df_tfbs_maf_es_dup=df_tfbs_maf_es_dup.join(df_tfbs_maf_es_dup_annotation.set_index(["Chromosome",
                                                                                  "Start",
                                                                                  "End",
                                                                                  "ID",
                                                                                  "REF",
                                                                                  "ALT",
                                                                                  "AF"]),
                                           on=["Chromosome","Start","End","ID","REF","ALT","AF"],
                                           how="left",
                                           rsuffix="_annotation")
df_tfbs_maf_es_dup.drop(columns=["promoter","enhancer","hepg2","k562"],inplace=True)
df_tfbs_maf_es_dup.rename(columns={"promoter_annotation":"promoter","enhancer_annotation":"enhancer","hepg2_annotation":"hepg2","k562_annotation":"k562"},inplace=True)
df_tfbs_maf_es_dup.drop_duplicates(subset=["Chromosome","Start","End","ID","REF","ALT","AF"],inplace=True)
df_tfbs_maf_es=pd.concat([df_tfbs_maf_es_nodup,df_tfbs_maf_es_dup],axis=0).reset_index(drop=True)

# determine transversion or transition
df_tfbs_maf_es["refalt"]=df_tfbs_maf_es["REF"]+df_tfbs_maf_es["ALT"]
transversions=["AG","GA","CT","TC"]
df_tfbs_maf_es["transversion"]=df_tfbs_maf_es["refalt"].isin(transversions).astype(int)
df_tfbs_maf_es.drop(columns=["refalt"],inplace=True)





#----------------------------
# train random forest 
#----------------------------

def train_rf(df,features,plot_suffix):
    x=df[features]
    y=df["Bin"]
    x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.2,random_state=42)
    rf=RandomForestClassifier(verbose=1,n_jobs=42)
    rf.fit(x_train,y_train)
    y_pred=rf.predict(x_test)
    # get feature importance
    feature_importance=pd.Series(rf.feature_importances_,index=features).sort_values(ascending=False)
    # plot feature importance
    plt.figure(figsize=(10,5))
    sns.barplot(x=feature_importance,y=feature_importance.index)
    # annotate with accuracy
    plt.text(0.75,0.05,f"accuracy: {(y_pred==y_test).mean()}",transform=plt.gca().transAxes)
    # annotate with MSE
    # plt.text(0.75,0.05,f"MSE: {((y_pred-y_test)**2).mean()}",transform=plt.g´ca().transAxes)
    plt.title("Feature importance")
    plt.savefig(f"rf_feature_importance_{plot_suffix}.pdf")
    plt.close()


# all
sub_tfs=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/sub_tfs.txt",header=None)[0].to_list()
super_tfs=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/super_tfs.txt",header=None)[0].to_list()
df_all=add_tf_annotation(df_tfbs_maf_es,sub_tfs,super_tfs)
# with sub/sper
features=[f"track_{i}" for i in range(8)]
features+=["score","promoter","enhancer","hepg2","k562",
          "transversion","sub-additive","super-additive","gini"]
train_rf(df_all,features,"all")

features=[f"track_{i}" for i in range(8)]
features+=["score","promoter","enhancer","hepg2","k562",
          "transversion","sub-additive","super-additive"]
train_rf(df_all,features,"all_no_gini")

# without sub/super
features=[f"track_{i}" for i in range(8)]
features+=["score","promoter","enhancer","hepg2","k562","transversion","gini"]
train_rf(df_all,features,"all_no_sub_super")

features=[f"track_{i}" for i in range(8)]
features+=["promoter","enhancer","hepg2","k562","transversion","sub-additive","super-additive"]
train_rf(df_all,features,"all_no_gini_score")



for cell_type in ["hepg2","k562"]:
    if cell_type=="hepg2":
        track_list=[0,2,4,6]
    if cell_type=="k562":
        track_list=[1,3,5,7]
    logger.info(f"train random forest for {cell_type}")
    df_subset=df_tfbs_maf_es[(df_tfbs_maf_es[cell_type]==1)].reset_index(drop=True)
    sub_tfs,super_tfs=get_sub_super_tfs(f"{cell_type}")
    df_subset=add_tf_annotation(df_subset,sub_tfs,super_tfs)
    # features=[f"track_{i}" for i in track_list]
    # features+=["score","transversion","sub-additive","super-additive","gini","promoter","enhancer"]
    # train_rf(df_subset,features,f"{cell_type}")
    # logger.info(f"train random forest for {cell_type}_no_sub_super")
    # features=[f"track_{i}" for i in track_list]
    # features+=["score","transversion","gini","promoter","enhancer",]
    # train_rf(df_subset,features,f"{cell_type}_no_sub_super")
    features=[f"track_{i}" for i in track_list]
    features+=["score","transversion","sub-additive","super-additive","promoter","enhancer"]
    train_rf(df_subset,features,f"{cell_type}_no_gini")








