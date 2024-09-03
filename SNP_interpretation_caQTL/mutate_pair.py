import pandas as pd
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer
from seq_annotators import JasparAnnotator, ReMapAnnotator
from seq_ops import SeqExtractor
from tf_cooperativity import write_pair_mutation,read_cooperativity


#---------------------------------
# Preprocess1: annotate with motif
#---------------------------------
df=pd.read_csv("Pd1_caQTL_annotated.csv")
df=df[df.log_10_Benjamini_Hochberg_Qvalue<=np.log10(0.05)].reset_index(drop=True)
df=df[['seqnames', 'start', 'end',"effect_size_.Pi.","txType"]]




# select tfs as remap
tfs_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/jaspar_tfs_hepg2.txt",header=None)[0].tolist()
tfs_split=split_dimer(tfs_hepg2)
tfs_hepg2_expressed=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Curate_motif_annot/expressed_tf_list_hepg2.tsv",sep="\t",header=None)[0].tolist()
tfs=tfs_hepg2+tfs_hepg2_expressed+tfs_split
tfs=list(set(tfs))



jaspar_annotator=JasparAnnotator(jaspar_path="/binf-isilon/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg19.bb")
remap_annotator=ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg19_Hep-G2.bed")


def annotate_overlapping_motif(df):
    best_motifs=[]
    best_scores=[]
    chip_evidences=[]
    for i in range(df.shape[0]):
        if i % 100 == 0:
            logger.info(f"Processing {i}th row")
        region=df.iloc[i,0:3]
        motifs=jaspar_annotator.annotate(region,by="1bp")
        if motifs.shape[0]==0:
            best_motifs.append("None")
            best_scores.append(0)
            chip_evidences.append(0)
            continue
        motifs=remap_annotator.annotate(motifs,region)
        motifs_sorted = motifs.sort_values(by=['chip_evidence', 'score'], ascending=[False, False])
        motif_region=motifs_sorted.iloc[0,0:3].tolist()
        best_motif=motifs_sorted.iloc[0,5]
        best_score=motifs_sorted.iloc[0,3]
        chip_evidence=motifs_sorted.iloc[0,6]
        best_motifs.append(best_motif)
        best_scores.append(best_score)
        chip_evidences.append(chip_evidence)
    df["motif"]=best_motifs
    df["score"]=best_scores
    df["chip_evidence"]=chip_evidences
    return df


# df=annotate_overlapping_motif(df)
# df.to_csv("Pd2_caQTL_motif_annotated.csv",index=False)


# select tfs expressed in HepG2
# df=pd.read_csv("Pd2_caQTL_motif_annotated.csv")
# # select tfs
# tfs_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/jaspar_tfs_hepg2.txt",header=None)[0].tolist()
# tfs_split=split_dimer(tfs_hepg2)
# tfs_hepg2_expressed=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Curate_motif_annot/expressed_tf_list_hepg2.tsv",sep="\t",header=None)[0].tolist()
# tfs=tfs_hepg2+tfs_hepg2_expressed+tfs_split
# tfs=list(set(tfs))
# # select rows with motif in tfs
# df=df[df["motif"].isin(tfs)].reset_index(drop=True)
# df.to_csv("Pd3_caQTL_motif_with_chip.csv",index=False)



#---------------------------------------------
# Preprocess2: annotate with additivity
#---------------------------------------------


# df=pd.read_csv("Pd3_caQTL_motif_with_chip.csv")
# df["start"]=df["end"]-300
# df["end"]=df["start"]+599


# # select tfs
# tfs_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/jaspar_tfs_hepg2.txt",header=None)[0].tolist()
# tfs_split=split_dimer(tfs_hepg2)
# tfs_hepg2_expressed=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Curate_motif_annot/expressed_tf_list_hepg2.tsv",sep="\t",header=None)[0].tolist()
# tfs=tfs_hepg2+tfs_hepg2_expressed+tfs_split
# tfs=list(set(tfs))


# # get tool
# seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg19.fa")
# jaspar_annotator=JasparAnnotator(jaspar_path="/binf-isilon/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg19.bb")
# remap_annotator=ReMapAnnotator(tfs)
# write_pair_mutation(df,seq_extractor,jaspar_annotator,remap_annotator,4,f"Pd4_caQTL_mutate_pair.csv")
# logger.info("Done")



#------------------------------------------------------------------------------
# Analysis: for same TF, compare its QTL in redundant v.s. codependent state
#------------------------------------------------------------------

df_motif=pd.read_csv("Pd3_caQTL_motif_with_chip.csv")
df_motif["region_idx"]="Region"+df_motif.index.astype(str)
# rename column "score" to "motif_score"
df_motif.rename(columns={"score":"motif_score"},inplace=True)
df_additivity=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/SNP_interpretation_caQTL/Pd4_caQTL_mutate_pair.csv",
                                 track_nums=[0,2,4,6])
# merge by region_idx
df=pd.merge(df_motif,df_additivity,on="region_idx",how="inner")

# select TFs with motif==protein2
df=df[df["motif"]==df["protein2"]].reset_index(drop=True)
df=df[df["motif_score"]==df["score2"]].reset_index(drop=True)
df=df[df["start_rel2"]<301].reset_index(drop=True)
df=df[df["end_rel2"]>300].reset_index(drop=True)




# group by 'seqnames','start','end','effect_size_.Pi.','txType','motif','motif_score', sum ci
df=df.groupby(['seqnames','start','end','effect_size_.Pi.','txType','motif','motif_score']).agg({"ci":"sum"}).reset_index()
df["codependency_state"]=(df["ci"]>0).astype(int)

df_sub=df[df["motif"]=="THRB"]
df_sub
# calculate pearsonr
from scipy.stats import pearsonr
pearsonr(df_sub["effect_size_.Pi."],df_sub["codependency_state"])

# group by motif and codependent state, calculate mean effect size
df=df.groupby(["motif","codependency_state"]).agg({"effect_size_.Pi.":"mean"}).reset_index()
# pivot table
df_pivot=df.pivot_table(index="motif",columns="codependency_state",values="effect_size_.Pi.").reset_index()
# remove rows with Nan
df_pivot=df_pivot.dropna().reset_index(drop=True)
# how many rows have column "0" > column "1"
df_pivot["diff"]=df_pivot[0]-df_pivot[1]
sum(df_pivot["diff"]<0)



# nohup python3 mutate_pair.py > mutate_pair.out &



#-------------------
# Archive
#-------------------




# df_qtl_promoter=df[df["txType"]=="promoter"].copy().reset_index(drop=True)
# df_qtl_enhancer=df[df["txType"]=="proximal"].copy().reset_index(drop=True)
# df_qtl_promoter=get_overlapping_motif(df_qtl_promoter)
# df_qtl_enhancer=get_overlapping_motif(df_qtl_enhancer)




# df_addtivity=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_additivity_property.csv")
# df_addtivity=df_addtivity[(df_addtivity["enhancers_hepg2"]=="super") & (df_addtivity["promoters_hepg2"]=="sub")].reset_index(drop=True)

# tfs=df_addtivity["protein"].unique()


# es_promoters=df_qtl_promoter[df_qtl_promoter["motif"].isin(tfs)].copy().reset_index(drop=True)
# es_promoters=es_promoters.groupby(["motif"]).agg({"effect_size_.Pi.":"mean"}).reset_index()
# es_promoters["re"]="promoter"
# es_enhancers=df_qtl_enhancer[df_qtl_enhancer["motif"].isin(tfs)].copy().reset_index(drop=True)
# es_enhancers=es_enhancers.groupby(["motif"]).agg({"effect_size_.Pi.":"mean"}).reset_index()
# es_enhancers["re"]="enhancer"

# # concat es_promoters and es_enhancers
# es=pd.concat([es_promoters,es_enhancers],ignore_index=True)
# # pivot table
# es_pivot=es.pivot_table(index="motif",columns="re",values="effect_size_.Pi.",aggfunc="mean").reset_index()
# # remove rows with Nan
# es_pivot=es_pivot.dropna().reset_index(drop=True)

# # wilcoxon signed
# from scipy.stats import wilcoxon
# wilcoxon(es_pivot["promoter"],es_pivot["enhancer"])