import pandas as pd


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_ratio, get_ppi



# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def read_all_files():
    # read cooperativity
    df_promoter_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_k562.csv",track_nums=[1,3,5,7])
    df_enhancer_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_k562.csv",track_nums=[1,3,5,7])
    df=pd.concat([df_promoter_k562,df_enhancer_k562],axis=0)
    # df=df_enhancer_k562
    return df




# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------
df=read_all_files()
df=calculate_tf_pair_cooperativity_ratio(df)
df_tf=df.groupby("protein2").agg({"redundancy_ci":"sum","codependency_ci":"sum"}).reset_index()
df_tf["sum_ci"]=df_tf["redundancy_ci"].abs()+df_tf["codependency_ci"]
df_tf=df_tf[df_tf["sum_ci"]>5].reset_index(drop=True)
df_tf["cooperativity_ratio"]=df_tf["codependency_ci"]/(df_tf["redundancy_ci"].abs()+df_tf["codependency_ci"])
df_tf.to_csv("tf_cooperativity_ratio_lenient.csv",index=False)

df=df[df["sum_ci"]>1].reset_index(drop=True)
df.to_csv("tf_pair_cooperativity_ratio_lenient.csv",index=False)


df=pd.read_csv("tf_pair_cooperativity_ratio_lenient.csv")
ppi_redundant=get_ppi(df,"redundancy",0.1,"redundancy_ppi.csv")
ppi_codependent=get_ppi(df,"codependency",0.9,"codependency_ppi.csv")


# df=pd.read_csv("tf_pair_cooperativity_ratio_enhancer_alone.csv")
# get_ppi(df,"redundancy",0.1,"redundancy_ppi_enhancer_alone.csv")
# get_ppi(df,"codependency",0.9,"codependency_ppi_enhancer_alone.csv")



# ----------------------------------------------------
# Inter v.s. intra community
# ----------------------------------------------------

def annotate_community(ppi,coop_graph,df_truth):
    """
    df_truth is the truth file for TF TF interaction
    """
    # create a dictionary of TFs and their community partners
    tf_community=dict()
    for community in ppi:
        for tf in community:
            if tf not in tf_community:
                tf_community[tf]=set()
            tf_community[tf]=tf_community[tf].union(community)
            tf_community[tf].remove(tf)
    # for each protein2 df_truth, check if it is in the same community as protein1
    res_intra=[]
    res_inter=[]
    for i in range(len(df_truth)):
        tf1=inweb["protein1"][i]
        tf2=inweb["protein2"][i]
        if tf2 in tf_community.keys() and tf1 in tf_community[tf2]:
            res_intra.append(True)
            res_inter.append(False)
        elif tf1 in coop_graph.index and tf2 in coop_graph.columns and coop_graph[tf1][tf2]==1:
            res_intra.append(False)
            res_inter.append(True)
        else:
            res_intra.append(False)
            res_inter.append(False)
    return res_intra,res_inter


inweb=pd.read_csv("/isdata/alab/people/pcr980/Resource/InWeb_TFCOF_PPIs.csv")
inweb=pd.read_csv("/isdata/alab/people/pcr980/Resource/STRING/String_TFCOF_PPIs.csv")
inweb=inweb[inweb["isProtein1TF"]].reset_index(drop=True)

df=pd.read_csv("tf_pair_cooperativity_ratio.csv")
df=df.pivot(index="protein1",columns="protein2",values="cooperativity_ratio")
res_intra,res_inter=annotate_community(ppi_codependent,(df>0.9).astype(int),inweb)
inweb["intra_codependent_community"]=res_intra
inweb["inter_codependent_community"]=res_inter

res_intra,res_inter=annotate_community(ppi_redundant,(df<0.1).astype(int),inweb)
inweb["intra_redundant_community"]=res_intra
inweb["inter_redundant_community"]=res_inter

# permutation of inweb for a null distribution
inweb_perm=inweb.copy()
inweb_perm["protein2"]=inweb_perm["protein2"].sample(frac=1).values
inweb_perm["same_codependent_community"]=annotate_community(ppi_codependent,inweb_perm)
inweb_perm["same_redundant_community"]=annotate_community(ppi_redundant,inweb_perm)


#--------------------------------------------------------------
# do codependent or redundant TFs interact with mediator more
#--------------------------------------------------------------
df_tf=pd.read_csv("tf_cooperativity_ratio_lenient.csv")
df_med=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/PPI_Mediator/2024-08-07_MED-TF_interactions.txt",sep="\t")
# how many genes in df_med are listed in df_tf
df_med["gene"].apply(lambda x: x in df_tf["protein2"].values).sum()
# df_med=df_med[df_med["significant"]==True].reset_index(drop=True)
# df_med["gene"].apply(lambda x: x in df_tf["protein2"].values).sum()

# select df with cooperativity_ratio>0.75
codependent_tfs=df_tf[df_tf["cooperativity_ratio"]>0.75]["protein2"].values
redundant_tfs=df_tf[df_tf["cooperativity_ratio"]<0.25]["protein2"].values

# how many codependent and redundant TFs interact with mediator
df_med["codependent"]=df_med["gene"].apply(lambda x: x in codependent_tfs)
df_med["redundant"]=df_med["gene"].apply(lambda x: x in redundant_tfs)
df_med["codependent"].sum()
df_med["redundant"].sum()

df_codependent=df_med[df_med["codependent"]]
df_redundant=df_med[df_med["redundant"]]

# fisher exact test
# is mediator more likely to interact with codependent or redundant TFs
from scipy.stats import fisher_exact
num_codependent_interactions=df_codependent["gene"].unique().shape[0]
num_codependent_noninteractions=len(codependent_tfs)-num_codependent_interactions
num_redundant_interactions=df_redundant["gene"].unique().shape[0]
num_redundant_noninteractions=len(redundant_tfs)-num_redundant_interactions
fisher_exact([[num_codependent_interactions,num_codependent_noninteractions],[num_redundant_interactions,num_redundant_noninteractions]])