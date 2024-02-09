import pandas as pd
import sys
import numpy as np

sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from region_ops import resize_df
from stat_tests import pearsonr_tolerating_nan
from seq_ops import SeqExtractor, resize_seq
from prediction import compute_predictions
# enhancer: 81046461-81046476
# promoter: 81107031-81107166
# not long, so always put mutated tile in the center


seq_extractor = SeqExtractor("/binf-isilon/alab/people/pcr980/Resource/hg19.fa")


#-------------------------------------
# clean .xlsx file, 
#-------------------------------------
df = pd.read_excel('media-1.xlsx', 
                   sheet_name='Supplementary Table 8',
                   skiprows=4)
# for columns Screen, grep "Enhancer" and "Promoter" and remove the rest content
df["Screen"]=df["Screen"].apply(lambda x: "Enhancer" if "Enhancer" in x else "Promoter" if "Promoter" in x else "None")
# split column 'VariantID' by ":", resulting in 3 columns named 'chromosome', 'start', 'mutation'
df[['chromosome','tile_start','mutation']]=df['VariantID'].str.split(":",expand=True)
# drop column 'VariantID'
df.drop(columns=['VariantID'], inplace=True)
# for column 'mutation', split by ">", resulting in 2 columns named 'REF' and 'ALT' 
df[['REF','ALT']]=df['mutation'].str.split(">",expand=True)
df["tile_end"]=df["tile_start"].astype(int)+df["REF"].str.len()
df["tile_start"]=df["tile_start"].astype(int)
df["tile_start"]=df["tile_start"]+1 # adjust for different coordinate system
df["tile_end"]=df["tile_end"]+1
assert np.all(df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["tile_start"], row["tile_end"]-1), axis=1)==df["REF"])

#-------------------------
# get reference sequence
#-------------------------
df['start']=df['tile_start']+(df['tile_end']-df['tile_start'])//2
df['end']=df['start']+1
df=resize_df(df,600,fix="center")
df["seq_ref"]=df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]-1), axis=1)


#-------------------------
# get alternative sequence
#-------------------------

# get prefix and suffix, make sure they are longer than 600bp in total 
df["prefix"]=df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["tile_start"]-350, row["tile_start"]-1), axis=1)
df["suffix"]=df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["tile_end"], row["tile_end"]+349), axis=1)

# assert df["prefix"] followed by df["REF"] followed by df["suffix"] contain df["seq_ref"] as a substring
assert np.all(df.apply(lambda row:  row["seq_ref"] in row["prefix"]+row["REF"]+row["suffix"], axis=1))

# change center region to ALT
df["seq_alt"]=df["prefix"]+df["ALT"]+df["suffix"]  
# resize to 600bp
df["seq_alt"]=df["seq_alt"].apply(lambda x: resize_seq(x)) 


# predict ref
pred_ref=compute_predictions(df["seq_ref"].tolist())

# predict alt
pred_alt=compute_predictions(df["seq_alt"].tolist())

# predicted delta
delta=pred_alt-pred_ref


#-------------------------
# calculate and plot correlation
#-------------------------

promoter_idx=(df["Screen"]=="Promoter")
enhancer_idx=(df["Screen"]=="Enhancer")

corrs=[]
pvals=[]
res=[]
tracks=[]
models=[]
for i in range(16):
    corr,pval=pearsonr_tolerating_nan(delta[promoter_idx,i],df["% change to PPIF expression"][promoter_idx])
    corrs.append(corr)
    pvals.append(pval)
    res.append(f"Promoter")
    tracks.append(i)
    if i>7:
        models.append("DeepCompare classification")
    else:
        models.append("DeepCompare regression")
    
for i in range(16):
    corr,pval=pearsonr_tolerating_nan(delta[enhancer_idx,i],df["% change to PPIF expression"][enhancer_idx])
    corrs.append(corr)
    pvals.append(pval)
    res.append(f"Enhancer")
    tracks.append(i)
    if i>7:
        models.append("DeepCompare classification")
    else:
        models.append("DeepCompare regression")
    
# get informatin of enformer:
for idx in [promoter_idx,enhancer_idx]:
    for col in ["Enformer CAGE TSS","Enformer DNase TSS ", "Enformer CAGE","Enformer DNase"]:
        corr,pval=pearsonr_tolerating_nan(df[col][idx],df["% change to PPIF expression"][idx])
        corrs.append(corr)
        pvals.append(pval)
        res.append(f"Promoter" if idx is promoter_idx else "Enhancer")
        tracks.append(-1)
        models.append(col)



df_corr=pd.DataFrame({"track":tracks,"correlation":corrs,"p-value":pvals,"re":res,"model":models})
df_corr.to_csv("correlation.csv",index=False)