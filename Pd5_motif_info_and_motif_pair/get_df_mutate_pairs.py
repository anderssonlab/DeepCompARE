import pandas as pd
import numpy as np
import pyranges as pr
import argparse
import sys
import os


sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_remap
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator, ReMapAnnotator
from loguru import logger
from tf_cooperativity import write_pair_mutation
#-----------------------
# Functions
#-----------------------


def analysis(file_prefix,device):
    remap_cell_type = get_remap(file_prefix)
    # load data and tools
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb")
    # method1: strict filtering of TFBS: have ChIP support
    # remap_annotator = ReMapAnnotator(f"/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_{remap_cell_type}.bed")
    # method2: loose filtering of TFBS: The TF is expressed
    tfs_expressed=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Curate_motif_annot/expressed_tf_list_hepg2.tsv",sep="\t",header=None).iloc[:,0].tolist()
    remap_annotator = ReMapAnnotator(tfs_expressed)
    
    # read the regions
    # TODO: change sep
    df_regions = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_prefix}.bed",sep="\t",header=None)
    # TODO: choose to -1
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    out_path = f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_lenient_{file_prefix}.csv"

    write_pair_mutation(df_regions,seq_extractor,jaspar_annotator,remap_annotator,device,out_path)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_prefix",type=str)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    
    analysis(args.file_prefix,args.device)
    logger.info(f"Done with {args.file_prefix}")