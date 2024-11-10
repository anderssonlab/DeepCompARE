import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from in_silico_mutagenesis import get_motif_isa
from gradxinp import compute_gradxinp
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator

#----------------
# Helper functions and tools
#----------------

seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")

jaspar_hepg2_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",score_thresh=500,
                                       rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv")

jaspar_k562_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",score_thresh=500,
                                      rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv")


def sample_regions(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_suffix}.bed",sep="\t",header=None).iloc[:,:3]
    df.columns=['chromosome','start','end']
    df=df.sample(500,axis=0).reset_index(drop=True)
    return df



def get_max_gradxinp(df_gradxinp,row,track_num):
    # Extract the max value in the specified range for the track
    gradxinps= df_gradxinp.loc[f"{row.seq_idx}_Track{track_num}", row.start_rel:row.end_rel].values
    return max(gradxinps)



def get_corr(df,jaspar_annotator,seq_extractor):
    motif_df=jaspar_annotator.annotate(df)
    motif_df=get_motif_isa(seq_extractor,motif_df)
    # calculate ism for each sequence
    df["sequences"]=df.apply(lambda row: seq_extractor.get_seq((row['chromosome'],row['start'],row['end']-1)),axis=1)
    df_gradxinp=compute_gradxinp(df["sequences"].tolist())
    df.drop(columns=['sequences'],inplace=True)
    # calculate 
    for track_num in range(16):
        motif_df[f'max_gradxinp_track{track_num}'] = motif_df.apply(lambda row: get_max_gradxinp(df_gradxinp, row, track_num), axis=1)
    corr_list=[]
    for track_num in range(16):
        corr,_=pearsonr(motif_df[f'isa_track{track_num}'],motif_df[f'max_gradxinp_track{track_num}'])
        corr_list.append(corr)
    return corr_list


#----------------
# Data generation
#----------------

corr_list=[]

for i in range(100):
    logger.info(f"iteration {i}")
    df_promoters_hepg2=sample_regions("promoters_hepg2")
    df_promoters_k562=sample_regions("promoters_k562")
    df_enhancers_hepg2=sample_regions("enhancers_hepg2")
    df_enhancers_k562=sample_regions("enhancers_k562")
    #
    df_hepg2=pd.concat([df_promoters_hepg2,df_enhancers_hepg2],ignore_index=True)
    df_k562=pd.concat([df_promoters_k562,df_enhancers_k562],ignore_index=True)
    #
    corr_hepg2=get_corr(df_hepg2,jaspar_hepg2_annotator,seq_extractor)
    corr_k562=get_corr(df_k562,jaspar_k562_annotator,seq_extractor)
    #
    corr_list.append(corr_hepg2)
    corr_list.append(corr_k562)


corr_df=pd.DataFrame(corr_list,columns=[f"track_{i}" for i in range(16)])
corr_df.to_csv("corr_motif_isa_vs_max_gradxinp.csv",index=False)


#---------------
# Plot
#---------------


# df=pd.read_csv("corr_motif_isa_vs_max_gradxinp.csv")
# df=df.melt(var_name="track",value_name="corr")
# # kde plot, hue=track
# plt.figure(figsize=(10,6))
# sns.kdeplot(data=df,x="corr",hue="track")
# plt.xlabel("Pearson correlation between motif ISA and max(gradxinp)")
# plt.savefig("corr_motif_isa_vs_max_gradxinp.pdf")
# plt.close()






# nohup python3 motif_isa_vs_max_gradxinp.py > motif_isa_vs_max_gradxinp.log &