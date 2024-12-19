#!/usr/bin/env python3
import pandas as pd
import torch
import os
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_annotators import JasparAnnotator

#-----------------------
# Functions
#-----------------------

THRESHOLD=500

def main(file_name):
    # Load data and tools
    df_regions=pd.read_csv(f"Pd1_Regions/{file_name}",sep="\t").iloc[:,0:3]
    df_regions.columns=["chromosome","start","end"]
    
    if "hepg2" in file_name:
        chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed"
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
    if "k562" in file_name:
        chip_file= "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed"
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                        score_thresh=THRESHOLD,
                                        chip_file=chip_file,
                                        rna_file=rna_file)

    jaspar_annotator.annotate(df_regions,outpath=f"Pd2_motif_info/motif_info_thresh_{THRESHOLD}_{file_name}.csv")

#-----------------------
# Analysis
#-----------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    
    args=parser.parse_args()
    main(args.file_name)
    logger.info(f"Done with {args.file_name}.")

