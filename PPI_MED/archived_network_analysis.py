import pandas as pd
import numpy as np

import seaborn as sns
from loguru import logger
import matplotlib.pyplot as plt
from collections import defaultdict

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_index, calculate_tf_cooperativity_index, make_symmetric, get_ppi


# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------
def process(cell_type):
    # read cooperativity
    df_promoter=read_cooperativity(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_{cell_type}.csv",track_nums=[1,3,5,7])
    df_enhancer=read_cooperativity(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_{cell_type}.csv",track_nums=[1,3,5,7])
    df=pd.concat([df_promoter,df_enhancer],axis=0)
    df=calculate_tf_pair_cooperativity_index(df)
    df.to_csv(f"tf_pair_cooperativity_index_{cell_type}_pre_filter.csv",index=False)
    # caution when filtering
    df=df[df["sum_ci"]>1].reset_index(drop=True)
    df.to_csv(f"tf_pair_cooperativity_index_{cell_type}_post_filter.csv",index=False)
    
    df_tf=calculate_tf_cooperativity_index(df)
    df_tf.to_csv(f"tf_cooperativity_index_{cell_type}.csv",index=False)





process("k562")
process("hepg2")


# ----------------------------------------------------
# 2. Get PPI
# ----------------------------------------------------
df=pd.read_csv("tf_pair_cooperativity_index_post_filter.csv")
ppi_redundant=get_ppi(df,"redundancy",0.1,"redundancy_ppi.csv")
ppi_codependent=get_ppi(df,"codependency",0.9,"codependency_ppi.csv")



# ----------------------------------------------------
# Inter v.s. intra community
# ----------------------------------------------------

def create_tf_community_members_dict(ppi):
    """
    ppi: a list of communities
    return a dictionary of TFs and their community partners
    """
    tf_community_members=dict()
    for community in ppi:
        for tf in community:
            if tf not in tf_community_members:
                tf_community_members[tf]=set()
            tf_community_members[tf]=tf_community_members[tf].union(community)
            tf_community_members[tf].remove(tf)
    return tf_community_members


def count_inter_intra_community_edges(ppi,coop_graph):
    tf_community_dict=create_tf_community_members_dict(ppi)
    all_nodes_ppi=set([item for sublist in ppi for item in sublist])
    edges_intra=0
    for tf,member_set in tf_community_dict.items():
        for member in member_set:
            edges_intra+=coop_graph[tf][member]
    edges_inter=0
    for tf,member_set in tf_community_dict.items():
        non_members=all_nodes_ppi.difference(member_set)
        for non_member in non_members:
            edges_inter+=coop_graph[tf][non_member]
    return edges_intra,edges_inter







def annotate_community(ppi,coop_graph,df_truth):
    """
    df_truth is the truth file for TF TF interaction
    """
    # create a dictionary of TFs and their community partners
    tf_community=create_tf_community_members_dict(ppi)
    # for each protein2 df_truth, check if it is in the same community as protein1
    res_intra=[]
    res_inter=[]
    for i in range(len(df_truth)):
        tf1=ppi_database["protein1"][i]
        tf2=ppi_database["protein2"][i]
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
string=pd.read_csv("/isdata/alab/people/pcr980/Resource/STRING/String_TFCOF_PPIs.csv")
ppi_database=pd.concat([inweb,string],axis=0).reset_index(drop=True)
ppi_database=ppi_database[ppi_database["isProtein1TF"]==True].reset_index(drop=True)

df=pd.read_csv("tf_pair_cooperativity_index_post_filter.csv")
df=make_symmetric(df)
ppi_database=ppi_database.merge(df[["protein1","protein2","cooperativity_index"]],on=["protein1","protein2"],how="left")


# annotate inter/intra community
df=df.pivot(index="protein1",columns="protein2",values="cooperativity_index")
codependency_graph=(df>0.9).astype(int)
res_intra,res_inter=annotate_community(ppi_codependent,codependency_graph,ppi_database)
ppi_database["intra_codependent_community"]=res_intra
ppi_database["inter_codependent_community"]=res_inter

redundancy_graph=(df<0.1).astype(int)
res_intra,res_inter=annotate_community(ppi_redundant,redundancy_graph,ppi_database)
ppi_database["intra_redundant_community"]=res_intra
ppi_database["inter_redundant_community"]=res_inter





# how many edges are in intra community/inter community?
ppi_database.intra_codependent_community.sum() # 501/1216
ppi_database.inter_codependent_community.sum() # 169/1215

ppi_database.intra_redundant_community.sum() # 575/1136
ppi_database.inter_redundant_community.sum() # 30/394


# permutation of ppi_database for a null distribution
ppi_database_perm=ppi_database.copy()
ppi_database_perm["protein2"]=ppi_database_perm["protein2"].sample(frac=1).values
ppi_database_perm["same_codependent_community"]=annotate_community(ppi_codependent,ppi_database_perm)
ppi_database_perm["same_redundant_community"]=annotate_community(ppi_redundant,ppi_database_perm)



count_inter_intra_community_edges(ppi_codependent,codependency_graph)
count_inter_intra_community_edges(ppi_redundant,redundancy_graph)



len(ppi_codependent)
len(ppi_redundant)


#--------------------------------------------------------------
# do codependent or redundant TF pairs interact with mediator more
#--------------------------------------------------------------
df_tf_pair=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_pair_cooperativity_index_post_filter.csv")
df_med=pd.read_csv("2024-08-07_MED-TF_interactions.txt",sep="\t")
# how many genes in df_med are listed in df_tf
df_med["gene"].apply(lambda x: x in df_tf["protein2"].values).sum()
df_med=df_med[df_med["significant"]==True].reset_index(drop=True)


df_tf_pair["protein1_interact_with_mediator"]=df_tf_pair["protein1"].apply(lambda x: x in df_med["gene"].values)
df_tf_pair["protein2_interact_with_mediator"]=df_tf_pair["protein2"].apply(lambda x: x in df_med["gene"].values)
df_tf_pair["both_interact_with_mediator"]=df_tf_pair["protein1_interact_with_mediator"] & df_tf_pair["protein2_interact_with_mediator"]
df_tf_pair["either_interact_with_mediator"]=df_tf_pair["protein1_interact_with_mediator"] | df_tf_pair["protein2_interact_with_mediator"]
df_tf_pair["neither_interact_with_mediator"]=~df_tf_pair["either_interact_with_mediator"]
# determine "both", "one", "neither" interact with mediator
df_tf_pair["interaction_with_mediator"]="neither"
df_tf_pair.loc[df_tf_pair["either_interact_with_mediator"],"interaction_with_mediator"]="one"
df_tf_pair.loc[df_tf_pair["both_interact_with_mediator"],"interaction_with_mediator"]="both"



# determine whether both interact with same subunit

# for each gene in df_med, create a set of baits

tf_bait_dict=defaultdict(set)
for row in df_med.itertuples():
    tf_bait_dict[row.gene].add(row.bait)
    
    

# df_tf_pair["same_subunit"]=df_tf_pair.apply(lambda x: len(tf_bait_dict[x["protein1"]].intersection(tf_bait_dict[x["protein2"]]))>0,axis=1)
# for df_tf_pair["same_subunit"]==True, set df_tf_pair["interaction_with_mediator"]="same_subunit"
# df_tf_pair.loc[df_tf_pair["same_subunit"],"interaction_with_mediator"]="same_subunit"

# plot distribution of cooperativity_index
sns.kdeplot(data=df_tf_pair,x="cooperativity_index",hue="interaction_with_mediator",common_norm=False)
plt.title("Cooperativity index distribution")
plt.savefig("cooperativity_index_interaction_with_mediator.pdf")
plt.close()

# mannwhitneyu test
from scipy.stats import mannwhitneyu
mannwhitneyu(df_tf_pair[df_tf_pair["interaction_with_mediator"]=="both"]["cooperativity_index"],df_tf_pair[df_tf_pair["interaction_with_mediator"]=="neither"]["cooperativity_index"])
mannwhitneyu(df_tf_pair[df_tf_pair["interaction_with_mediator"]=="one"]["cooperativity_index"],df_tf_pair[df_tf_pair["interaction_with_mediator"]=="neither"]["cooperativity_index"])
mannwhitneyu(df_tf_pair[df_tf_pair["interaction_with_mediator"]=="both"]["cooperativity_index"],df_tf_pair[df_tf_pair["interaction_with_mediator"]=="one"]["cooperativity_index"])
#mannwhitneyu(df_tf_pair[df_tf_pair["interaction_with_mediator"]=="same_subunit"]["cooperativity_index"],df_tf_pair[df_tf_pair["interaction_with_mediator"]=="both"]["cooperativity_index"])

df_tf_pair[df_tf_pair["interaction_with_mediator"]=="same_subunit"]["cooperativity_index"].describe()
df_tf_pair[df_tf_pair["interaction_with_mediator"]=="both"]["cooperativity_index"].describe()

