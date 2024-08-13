import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator, ReMapAnnotator
from tf_cooperativity import write_pair_mutation, read_cooperativity, calculate_tf_pair_cooperativity_ratio, get_ppi


#-----------------------
# get df_mutate_pair
#-----------------------
def get_tfs():
    df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/PPI_K562/2024-08-07_MED-TF_interactions.txt",sep="\t")
    # select df["significant"] == true
    df=df[df["significant"]==True].reset_index(drop=True)
    return df["gene"].unique().tolist()


def get_df_mutate_pair(out_path,device=7):
    # load data and tools
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb")
    remap_annotator = ReMapAnnotator(get_tfs())
    
    # read the regions
    df_promoters= pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/promoters_k562.bed",sep="\t",header=None)
    df_enhancers= pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/enhancers_k562.bed",sep="\t",header=None)
    df_regions=pd.concat([df_promoters,df_enhancers],axis=0)
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    write_pair_mutation(df_regions,seq_extractor,jaspar_annotator,remap_annotator,device,out_path)


# get_df_mutate_pair("df_mutate_pairs_MED_tfs.csv")
# logger.info(f"Done!")


# nohup python3 get_df_mutate_pairs.py > get_df_mutate_pairs.log 2>&1 & 


#-----------------------
# Analysis
#-----------------------


df=read_cooperativity("df_mutate_pairs_MED_tfs.csv",track_nums=[1,3,5,7])
df=calculate_tf_pair_cooperativity_ratio(df)
df.to_csv("tf_pair_cooperativity_ratio.csv",index=False)
get_ppi(df,"redundancy",0.1,"redundancy_ppi_MED_tfs.csv")
get_ppi(df,"codependency",0.95,"codependency_ppi_MED_tfs.csv")




df=pd.read_csv("tf_pair_cooperativity_ratio.csv")
# remove rows where protein1==protein2
df=df[df["protein1"]!=df["protein2"]].reset_index(drop=True)
df=df[df["sum_ci"]>0.5].reset_index(drop=True)
df_med=pd.read_csv("2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]==True].reset_index(drop=True)
# creat a dict: key is gene, value is bait
gene_bait_dict={}
for tf in df_med.gene.unique():
    gene_bait_dict[tf]=df_med[df_med["gene"]==tf]["bait"].unique().tolist()
    
# for each row in df
# count baits shared by protein1 and protein2 

def count_shared_baits(row):
    baits1=gene_bait_dict[row["protein1"]]
    baits2=gene_bait_dict[row["protein2"]]
    return len(set(baits1).intersection(set(baits2)))

df["share_baits"]=df.apply(count_shared_baits,axis=1)
# df["share_baits"]=(df["share_baits"]>0).astype(int)


# fisher exact test: are cooperativity ratio <0.3 enriched in share_baits==1?
from scipy.stats import fisher_exact
df["cooperativity_ratio<0.2"]=df["cooperativity_ratio"]<0.2
df["cooperativity_ratio<0.2"]=df["cooperativity_ratio<0.2"].astype(int)
contingency_table=pd.crosstab(df["cooperativity_ratio<0.2"],df["share_baits"])
oddsratio, pvalue = fisher_exact(contingency_table)



# plot distribution of cooperativity ratio, hue is share_baits
sns.kdeplot(x="cooperativity_ratio",hue="share_baits",data=df)
plt.savefig("share_baits_vs_cooperativity_ratio.png")
plt.close()