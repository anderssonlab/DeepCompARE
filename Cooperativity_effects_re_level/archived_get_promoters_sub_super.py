import pandas as pd



cell_type="hepg2"

df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_{cell_type}_track0.csv")
# select true chip
df=df[df["chip_evidence"]].reset_index(drop=True)
# add sub/super annotation

sub_tfs=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/sub_tfs.txt",header=None)[0].tolist()
super_tfs=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/super_tfs.txt",header=None)[0].tolist()


# annotate column "protein"
df["tf_type"]="unknown"
df.loc[df["protein"].isin(sub_tfs),"tf_type"]="sub"
df.loc[df["protein"].isin(super_tfs),"tf_type"]="super"


df_grouped=df.groupby("seq_idx")["tf_type"].value_counts().unstack().fillna(0)


# read promoter ranges
promoter_seqs=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/promoters_{cell_type}.bed",sep="\t",header=None)
# label promoters i with f"Seq{i}"
promoter_seqs["seq_idx"]=promoter_seqs.index    
promoter_seqs["seq_idx"]=promoter_seqs["seq_idx"].apply(lambda x: f"Seq{x}")

# merge with df_grouped
df_grouped=df_grouped.reset_index()
df_grouped=pd.merge(promoter_seqs,df_grouped,on="seq_idx",how="left")
df_grouped.to_csv(f"promoters_sub_super_counts_{cell_type}.csv",header=True,index=False)