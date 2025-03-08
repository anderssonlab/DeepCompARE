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
    df_regions=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/DCMain_Fig4g_cooperativity_tf_location/Pd1_E1E2P1P2/{file_name}.tsv",sep="\t").iloc[:,0:3]
    df_regions.columns=["chromosome","start","end"]
    # select chromosomes
    df_regions=df_regions[df_regions["chromosome"].apply(lambda x: x in [f"chr{i}" for i in list(range(1,23))+["X","Y"]])].reset_index(drop=True)
    if "k562" in file_name:
        chip_file= "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg19_k562.bed"
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg19.bb",
                                        score_thresh=THRESHOLD,
                                        chip_file=chip_file,
                                        rna_file=rna_file)
    jaspar_annotator.annotate(df_regions,outpath=f"motif_info_thresh_{THRESHOLD}_{file_name}.csv")



#-----------------------
# Analysis
#-----------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    
    args=parser.parse_args()
    main(args.file_name)
    logger.info(f"Done with annotating {args.file_name}.")

