import pandas as pd
import sys

sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import pearsonr_tolerating_nan
from seq_ops import SeqExtractor
from prediction import compute_predictions



seq_extractor = SeqExtractor("/binf-isilon/alab/people/pcr980/Resource/hg19.fa")

df = pd.read_excel('OddsRatioCBEABEdCas9.xlsx',index_col=0)
df_info=pd.DataFrame({"chromosome":"chr12",
                      "start":df.iloc[:,0],
                      "end":df.iloc[:,0]+1})
df_info["prefix"]=df_info.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"]-299, row["start"]-1), axis=1)
df_info["suffix"]=df_info.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["end"], row["end"]+299), axis=1)
df_info["ref_base"]=df_info.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]-1), axis=1)
df_info["seq_ref"]=df_info["prefix"]+df_info["ref_base"]+df_info["suffix"]
df_info["alt_base"]="N"
df_info["seq_alt"]=df_info["prefix"]+df_info["alt_base"]+df_info["suffix"]
df_pred_ref=compute_predictions(df_info["seq_ref"].tolist())
df_pred_alt=compute_predictions(df_info["seq_alt"].tolist())

delta=df_pred_alt-df_pred_ref
delta=pd.DataFrame(delta)
delta.index=df_info.index

colnames=[]
for model in ["reg","class"]:
    for modality in ["cage", "dhs","starr","sure"]:
        for cell in ["hepg2","k562"]:
            colnames.append(f"{model}_{modality}_{cell}")
            
delta.columns=colnames
delta.to_csv("delta.csv")
pearsonr_tolerating_nan(df.iloc[:,7],df.iloc[:,8])

for i in range(16):
    print(pearsonr_tolerating_nan(delta.iloc[:,i],df.iloc[:,7:9].mean(axis=1)))

