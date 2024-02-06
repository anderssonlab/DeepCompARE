import pandas as pd
import pyBigWig
import pyranges as pr
from loguru import logger  
import argparse


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import SeqExtractor
from region_ops import calc_gc_context



# function
def annotate_region(region,tf,chip_evidence):
    """
    For each region, get jaspar information for tf only.
    Return a data frame, with chip_evidence column added
    
    """
    try:
        jaspar_info_list = jaspar.entries(region[0], region[1], region[2])
    except:
        return pd.DataFrame(columns=['start', 'end', 'protein', 'score', 'strand', 'chromosome','chip_evidence', 'motif_seq', 'motif_gc', 'context_gc'])

    df = pd.DataFrame(jaspar_info_list, columns=['start', 'end', 'details'])
    if df.shape[0]==0:
        return pd.DataFrame(columns=['start', 'end', 'protein', 'score', 'strand', 'chromosome','chip_evidence', 'motif_seq', 'motif_gc', 'context_gc'])
    df[['protein', 'score', 'strand']] = df['details'].str.split('\t', expand=True)
    df.drop(columns='details', inplace=True)
    df.protein=df.protein.str.upper()
    df=df[df.protein==tf]
    if df.shape[0]==0:
        return pd.DataFrame(columns=['start', 'end', 'protein', 'score', 'strand', 'chromosome','chip_evidence', 'motif_seq', 'motif_gc', 'context_gc'])
    df["chromosome"]=region[0]
    df["chip_evidence"]=chip_evidence
    df["motif_seq"]=df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
    df["motif_gc"]=calc_gc_context(df,0,seq_extractor)
    df["context_gc_2bp"]=calc_gc_context(df,2,seq_extractor)
    df["context_gc_10bp"]=calc_gc_context(df,10,seq_extractor)
    df["context_gc_50bp"]=calc_gc_context(df,50,seq_extractor)
    df["context_gc_100bp"]=calc_gc_context(df,100,seq_extractor)
    df["context_gc_300bp"]=calc_gc_context(df,300,seq_extractor)
    return df


parser=argparse.ArgumentParser()
parser.add_argument("--tf",type=str,help="TF name")
tf=parser.parse_args().tf

# read in jaspar
jaspar_path="/binf-isilon/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb"
jaspar = pyBigWig.open(jaspar_path)
jaspar.isBigBed()


# read in open chromatin regions
gr_ocr=pr.read_bed("/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/DHS_K562_loose.bed")
gr_ocr=gr_ocr.merge()


# read in ChIP peak regions
df_chip=pd.read_csv("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed",sep="\t",header=None)
df_chip.columns=["chromosome","start","end","protein","score","strand","thick_start","thick_end","rgb"]
df_chip=df_chip[df_chip["chromosome"].isin([f"chr{i}" for i in range(1,23)]+["chrX","chrY"])].reset_index(drop=True)
df_chip.protein=df_chip.protein.str.split(":",expand=True)[0]
df_chip.protein=df_chip.protein.str.upper()
df_chip.columns=["Chromosome","Start","End","protein","score","strand","thick_start","thick_end","rgb"]
df_chip=df_chip[df_chip.protein==tf]

gr_chip_pos=pr.PyRanges(df_chip)
gr_chip_neg=gr_ocr.subtract(gr_chip_pos)

df_chip_pos=gr_chip_pos.as_df()
df_chip_neg=gr_chip_neg.as_df()


seq_extractor = SeqExtractor()

for i in range(df_chip_pos.shape[0]):
    if i%1000==0:
        logger.info(f"Processing {i}th row of positve set")
    region=df_chip_pos.iloc[i,0:3].tolist()
    df_region=annotate_region(region,tf,"True")
    if i==0:
        df_region.to_csv(f"Pd1_motif_loc_and_context_gc/{tf}.csv",index=False,mode="w")
    else:
        df_region.to_csv(f"Pd1_motif_loc_and_context_gc/{tf}.csv",index=False,mode="a",header=False)
    
for i in range(df_chip_neg.shape[0]):
    if i%1000==0:
        logger.info(f"Processing {i}th row of negative set")
    region=df_chip_neg.iloc[i,0:3].tolist()
    df_region=annotate_region(region,tf,"False")
    df_region.to_csv(f"Pd1_motif_loc_and_context_gc/{tf}.csv",index=False,mode="a",header=False)
        