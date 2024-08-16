import seaborn as sns
from loguru import logger
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr,pearsonr

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer



# TODO: remove proteins rarely studied (Dewei)


#----------------------------
# Analysis 1: Get score for all tfs
#----------------------------
def tf_analysis(network, protein, threshold):
    with open(network,'rt') as f:
        reader = csv.reader(f, delimiter=' ')
        next(reader)
        links = set()
        for p1,p2,score in reader:
            if protein in (p1,p2) and float(score)/1000 >= threshold:
                # sort the p1 and p2
                if p1 == protein:
                    links.add('\t'.join([p1,p2,str(float(score)/1000)]))
                else:
                    links.add('\t'.join([p2,p1,str(float(score)/1000)]))
        score_sum = sum([float(x.split('\t')[-1]) for x in links])
    return protein,score_sum


# get string id for all tfs
# all columns should be string
aliases = pd.read_csv('/isdata/alab/people/pcr980/Resource/STRING/9606.protein.aliases.v12.0.txt', sep='\t', dtype=str)

df_tf_pair=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_pair_cooperativity_ratio_pre_filter_lenient.csv")
all_tfs=set(list(df_tf_pair['protein2'].unique())+list(df_tf_pair['protein2'].unique()))
all_tfs=split_dimer(all_tfs)
df_mapper=aliases[aliases['alias'].isin(all_tfs)].iloc[:,:2].drop_duplicates().reset_index(drop=True)




# perform tf_analysis for all tfs
results_physical = list()
network='/isdata/alab/people/pcr980/Resource/STRING/9606.protein.physical.links.v12.0.txt'
for id in df_mapper['#string_protein_id']:
    print(id)
    score_sum = tf_analysis(network, id, 0.7)
    results_physical.append(score_sum)


results_all_links = list()
network='/isdata/alab/people/pcr980/Resource/STRING/9606.protein.links.v12.0.txt'
for id in df_mapper['#string_protein_id']:
    print(id)
    score_sum = tf_analysis(network, id, 0.7)
    results_all_links.append(score_sum)


# convert results to dataframe
df_physical=pd.DataFrame(results_physical,columns=['id','score_physical'])
df_all_links=pd.DataFrame(results_all_links,columns=['id','score_all_links'])
# merge by id
df_results=df_physical.merge(df_all_links,on='id',how='left')
# merge with df_mapper
df_results=df_results.merge(df_mapper,left_on='id',right_on='#string_protein_id',how='left')
df_results.to_csv('string_protein_interaction_scores_all_tfs.csv',index=False)



#--------------------------------------------
# Analysis 2. analyze redundant v.s. codependent interaction score
#--------------------------------------------
tfs_codependent = "/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_codependent_lenient.txt"
tfs_codependent = open(tfs_codependent).read().strip().split('\n')
tfs_redundant = "/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_redundant_lenient.txt"
tfs_redundant = open(tfs_redundant).read().strip().split('\n')

df_results=pd.read_csv('string_protein_interaction_scores_all_tfs.csv')
df_codependent=df_results[df_results['alias'].isin(tfs_codependent)].groupby('alias').agg({'score_physical':'sum','score_all_links':'sum'}).reset_index()
df_codependent['tf_type']='codependent'
df_redundant=df_results[df_results['alias'].isin(tfs_redundant)].groupby('alias').agg({'score_physical':'sum','score_all_links':'sum'}).reset_index()
df_redundant['tf_type']='redundant'





tfs=pd.concat([df_codependent,df_redundant])
sns.kdeplot(data=tfs,x='score_physical',hue='tf_type',common_norm=False)
plt.title("STRING Protein interaction score (physical interaction)")
plt.savefig('Plots/string_score_physical_lenient.pdf')
plt.close()

sns.kdeplot(data=tfs,x='score_all_links',hue='tf_type',common_norm=False)
plt.title("STRING Protein interaction score (all links)")
plt.savefig('Plots/string_score_all_links_lenient.pdf')
plt.close()

# mannwhitneyu test
mannwhitneyu(df_codependent['score_physical'],df_redundant['score_physical'])
mannwhitneyu(df_codependent['score_all_links'],df_redundant['score_all_links'])

df_codependent['score_physical'].median()
df_redundant['score_physical'].median()

df_codependent['score_all_links'].median()
df_redundant['score_all_links'].median()


#------------------------------------------------------------
# spearman correlation between cooperativity ratio and protein interaction
#------------------------------------------------------------
df_cooperativity_ratio=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_cooperativity_ratio_lenient.csv")
# merge with df_results
df_results=df_results.merge(df_cooperativity_ratio,left_on='alias',right_on='protein2',how='inner')
spearmanr(df_results['cooperativity_ratio'],df_results['score_physical'])
spearmanr(df_results['cooperativity_ratio'],df_results['score_all_links'])



#-----------------------
# TF pair analysis
#-----------------------

# create mapper
aliases = pd.read_csv('/isdata/alab/people/pcr980/Resource/STRING/9606.protein.aliases.v12.0.txt', sep='\t', dtype=str)
df_tf_pair=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_pair_cooperativity_ratio_post_filter_lenient.csv")
all_tfs=list(df_tf_pair['protein2'].unique())
all_tfs=split_dimer(all_tfs)
df_mapper=aliases[aliases['alias'].isin(all_tfs)].iloc[:,:2].drop_duplicates().reset_index(drop=True)
# are all_tfs in the aliases? Yes
len(all_tfs)==len(df_mapper["alias"].unique())
# are df_mapper["#string_protein_id"] unique? No
len(df_mapper["#string_protein_id"].unique())==len(df_mapper["#string_protein_id"])
# get rows where "#string_protein_id" is not unique 
# same string ID: NFE2L1-NRF1 HNF1A-HNF4a
df_mapper["#string_protein_id"].value_counts()
df_mapper[df_mapper["#string_protein_id"]=="9606.ENSP00000354855"]


df_network=pd.read_csv('/isdata/alab/people/pcr980/Resource/STRING/9606.protein.links.v12.0.txt',sep=' ')
# remove rows with protein1 not in df_mapper["#string_protein_id"]
df_network=df_network[df_network["protein1"].isin(df_mapper["#string_protein_id"])].reset_index(drop=True)
# remove rows with protein2 not in df_mapper["#string_protein_id"]
df_network=df_network[df_network["protein2"].isin(df_mapper["#string_protein_id"])].reset_index(drop=True)

# replace protein1 and protein2 with protein names
df_network=df_network.merge(df_mapper,left_on="protein1",right_on="#string_protein_id",how="left")
df_network=df_network.merge(df_mapper,left_on="protein2",right_on="#string_protein_id",how="left")
df_network=df_network.drop(["protein1","protein2","#string_protein_id_x","#string_protein_id_y"],axis=1)
# rename alias_x and alias_y
df_network.columns=["score","protein1","protein2"]

# merge with df_tf_pair by protein1 and protein2
df_network=df_network.merge(df_tf_pair,left_on=["protein1","protein2"],right_on=["protein1","protein2"],how="left")
# drop rows with NaN
df_network=df_network.dropna().reset_index(drop=True)
# group by protein1 and protein2, select max "cooperativity_ratio" and mean "score"
df_network=df_network.groupby(["protein1","protein2"]).agg({"cooperativity_ratio":"max","score":"mean"}).reset_index()


# spearman correlation between cooperativity_ratio and score
from scipy.stats import spearmanr
spearmanr(df_network["cooperativity_ratio"],df_network["score"]) # corr=0.08 p<1e-4

# scatter plot the cooperativity_ratio and score
import seaborn as sns
sns.scatterplot(data=df_network,x="cooperativity_ratio",y="score")
plt.xlabel("Cooperativity ratio")
plt.ylabel("Protein interaction score")
plt.title("Cooperativity ratio and protein interaction score")
plt.savefig("Plots/cooperativity_ratio_protein_interaction_score_lenient.pdf")
plt.close()