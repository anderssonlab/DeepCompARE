import pandas as pd



# get tfs that are redundant in promoter and codependent in enhancer
df_cooperativity=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_cooperativity_ratio_promoters_vs_enhancers.csv")
# select rows with cr_promoter<0.3 and cr_enhancer>0.7
df_target=df_cooperativity[(df_cooperativity['cr_promoters']<0.3) & (df_cooperativity['cr_enhancers']>0.7)]
df_target2=df_cooperativity[(df_cooperativity['cr_promoters']>0.7) & (df_cooperativity['cr_enhancers']<0.3)]






df_promoter_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_promoters_hepg2.csv")
df_promoter_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_promoters_k562.csv")
df_enhancer_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_enhancers_hepg2.csv")
df_enhancer_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_enhancers_k562.csv")
df_promoter_hepg2["dataset"]="promoter_hepg2"
df_promoter_k562["dataset"]="promoter_k562"
df_enhancer_hepg2["dataset"]="enhancer_hepg2"
df_enhancer_k562["dataset"]="enhancer_k562"
df=pd.concat([df_promoter_hepg2,df_promoter_k562,df_enhancer_hepg2,df_enhancer_k562])
df["distance"]=(df["start_rel2"]-df["start_rel1"]).abs()
# remove rows with distance larger than 200
df=df[df["distance"]<200].reset_index(drop=True)
df["re"]=df["dataset"].apply(lambda x: x.replace("_hepg2","").replace("_k562",""))





# group by dataset, protein1, protein2, and count the number of rows
df_grouped=df.groupby(["re","protein1","protein2"]).size().reset_index(name="count")
# select rows where protein2 is in df_target["protein2"]
df_grouped=df_grouped[df_grouped["protein2"].isin(df_target["protein2"])].reset_index(drop=True)
df_grouped=df_grouped.sort_values(by=["protein2","re"])

# for each protein2 and re, get top 10 protein1 sorted by count
df_grouped=df_grouped.groupby(["protein2","re"]).apply(lambda x: x.nlargest(10,"count")).reset_index(drop=True)
df_grouped.to_csv("top10_protein1_by_count.csv",index=False)