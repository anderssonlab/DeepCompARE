import pandas as pd
import torch
import os
import pyranges as pr
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from prediction import compute_predictions

#-----------------------
# Functions
#-----------------------



def main(file_name,device):
    # Load data and tools
    df=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_name}.bed").df.iloc[:,0:3]
    df.columns=["chromosome","start","end"]
    df["end"]=df["end"]-1
    
    # step 2: get sequence
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    seqs=df.apply(lambda row: seq_extractor.get_seq((row["chromosome"],row["start"],row["end"])),axis=1)
    preds=compute_predictions(seqs,device=device)
    pd.DataFrame(preds).to_csv(f"predictions_{file_name}.csv",index=False)


#-----------------------
# Analysis
#-----------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    parser.add_argument('--device', type=str)
    
    args=parser.parse_args()
    main(args.file_name,args.device)
    logger.info(f"Finished annotating {args.file_name}.")

