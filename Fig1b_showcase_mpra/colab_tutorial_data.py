import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import pearsonr
import numpy as np

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator



#------------------------------
# Load annotators
#------------------------------
seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
jaspar_hepg2_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                       chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed",
                                       rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
                                       )
jaspar_k562_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                       chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed",
                                       rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
                                       )


#------------------------------
# write target sequences: 
#------------------------------
# F9 promoter: chrX:139530463-139530765
# LDLR promoter: chr19:11,089,231-11,089,548
# PKLR promoter: chr1:155,301,395-155,301,864
# SORT1 enhancer: chr1:109274652-109275251


df = pd.DataFrame([
    {"cell_line":"HepG2","element_name": "F9 promoter", "chrom": "chrX", "start": 139530463, "end": 139530765},
    {"cell_line":"HepG2","element_name": "LDLR promoter", "chrom": "chr19", "start": 11089231, "end": 11089548},
    {"cell_line":"K562","element_name": "PKLR promoter", "chrom": "chr1", "start": 155301395, "end": 155301864},
    {"cell_line":"HepG2","element_name": "SORT1 enhancer", "chrom": "chr1", "start": 109274652, "end": 109275251},
])
df["sequence"] = df.apply(lambda row: seq_extractor.get_seq((row["chrom"], row["start"], row["end"])), axis=1)

df.to_csv("mpra_region_sequences.csv", index=False)



#-------------------------------
# get motif list for all 4 regions
#-------------------------------


def reduce_protein_names(protein_list):
    protein_list=list(set(protein_list))
    # order alphabetically
    protein_list.sort()
    # if there are more than 2 proteins sharing same prefix of length > 4
    # only keep the prefix, followed by "s"
    # eg: hoxa9, hoxa9b, hoxa9c -> hoxa9s
    protein_dict={}
    for protein in protein_list:
        prefix=protein[:4]
        if prefix in protein_dict:
            protein_dict[prefix].append(protein)
        else:
            protein_dict[prefix]=[protein]
    prefix_list=[]
    for prefix in protein_dict:
        if len(protein_dict[prefix])>1:
            prefix_list.append(prefix+"s")
        else:
            prefix_list.append(protein_dict[prefix][0])
    # concatenate by "\n"
    return "\n".join(prefix_list)






def reduce_motifs(df_motif,window=4):
    # if start is within 3bp of another start
    # choose the top 3 based on "score"
    # concatenate protein with "\n", use largest "end" as end
    df_res=pd.DataFrame()
    while df_motif.shape[0]>0:
        current_start=df_motif.loc[0,"start_rel"]
        df_temp=df_motif[(df_motif["start_rel"]>=current_start) & (df_motif["start_rel"]<=(current_start+window))].copy().reset_index(drop=True)
        df_temp=df_temp.sort_values(by="score",ascending=False).reset_index(drop=True)
        df_temp=df_temp.iloc[:2,:]
        df_temp["protein"]=reduce_protein_names(df_temp["protein"])
        df_temp["end_rel"]=df_temp["end_rel"].max()
        df_temp["start_rel"]=df_temp["start_rel"].min()
        df_temp["start"]=df_temp["start"].min()
        df_temp["end"]=df_temp["end"].max()
        df_res=df_res.append(df_temp.iloc[0,:],ignore_index=True)
        # remove the rows in df_temp from df_motif
        df_motif=df_motif[df_motif["start_rel"]>current_start+window].copy().reset_index(drop=True)
    return df_res




def get_motifs(jaspar_annotator,region,score_threshold):
    df_motif=jaspar_annotator.annotate(region)
    # remove uncertain proteins
    df_motif=df_motif[~df_motif["protein"].str.contains("::")].copy().reset_index(drop=True)
    df_chip=df_motif.loc[df_motif["chip_evidence"]==True,:].reset_index(drop=True)
    df_rest=df_motif.loc[df_motif["chip_evidence"]==False,:].reset_index(drop=True)
    df_rna=df_rest.loc[df_rest["rna_evidence"]==True,:].reset_index(drop=True)
    df_rna=df_rna.loc[df_rna["score"]>=score_threshold,:].reset_index(drop=True)
    df_motif=pd.concat([df_chip,df_rna],axis=0).reset_index(drop=True)
    # sort df_motif by "start"
    df_motif=df_motif.sort_values(by="start").reset_index(drop=True)
    df_motif["start_rel"]=df_motif["start"]-region[1]
    df_motif["end_rel"]=df_motif["end"]-region[1]
    df_motif=reduce_motifs(df_motif)
    return df_motif



def extract_motifs(row):
    region = (row["chrom"], row["start"], row["end"])
    cell_line = row["cell_line"]
    if cell_line == "HepG2":
        jaspar_annotator = jaspar_hepg2_annotator
    elif cell_line == "K562":
        jaspar_annotator = jaspar_k562_annotator
    else:
        raise ValueError(f"Unsupported cell line: {cell_line}")
    df_motif = get_motifs(jaspar_annotator, region, score_threshold=360)
    df_motif["element_name"] = row["element_name"]
    return df_motif









# Concatenate motif annotations
df_motifs = pd.concat(df.apply(extract_motifs, axis=1).to_list(), ignore_index=True)

# Save to CSV
df_motifs.to_csv("motif_annotations.csv", index=False)
