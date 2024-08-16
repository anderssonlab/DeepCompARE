import pandas as pd
import argparse
import sys


sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_binding_evidence
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator, BindingEvidenceAnnotator
from loguru import logger
from tf_cooperativity import write_pair_mutation
#-----------------------
# Functions
#-----------------------


def analysis(file_prefix,device,sep="\t"):
    be_type = get_binding_evidence(file_prefix)
    # load data and tools
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",by="contained")
    # method1: strict filtering of TFBS: have ChIP support
    # be_annotator = BindingEvidenceAnnotator(f"/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_{be_type}.bed")
    # method2: loose filtering of TFBS: The TF is expressed
    be_annotator = BindingEvidenceAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_{be_type}.tsv",mode="rna")
    
    df_regions = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_prefix}.bed",sep=sep,header=None)
    # TODO: choose to -1 or not
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    out_path = f"mutate_pairs_lenient_{file_prefix}.csv"
    write_pair_mutation(df_regions,seq_extractor,jaspar_annotator,be_annotator,device,out_path)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_prefix",type=str)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    analysis(args.file_prefix,args.device)
    logger.info(f"Done with {args.file_prefix}")