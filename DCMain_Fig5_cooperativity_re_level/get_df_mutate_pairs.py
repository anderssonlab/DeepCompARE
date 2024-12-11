import pandas as pd
import argparse
import sys


sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator
from loguru import logger
from tf_cooperativity import write_pair_mutation

#-----------------------
# Functions
#-----------------------

GENOME="hg38"

def analysis(file_name,device):
    # load data and tools
    seq_extractor = SeqExtractor(f"/isdata/alab/people/pcr980/Resource/{GENOME}.fa")
    if "hepg2" in file_name:
        jaspar_annotator=JasparAnnotator(f"/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_{GENOME}.bb",
                                        score_thresh=500,
                                        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv",
                                        chip_file=f"/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_{GENOME}_hepg2.bed",
                                        subset="rna")
    elif "k562" in file_name:
        jaspar_annotator=JasparAnnotator(f"/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_{GENOME}.bb",
                                        score_thresh=500,
                                        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv",
                                        chip_file=f"/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_{GENOME}_k562.bed",
                                        subset="rna")
    else:
        raise ValueError("file_name should contain hepg2 or k562")
    df_regions = pd.read_csv(f"Pd1_Regions/{file_name}.tsv",sep="\t")
    #df_regions = pd.read_csv(f"Pd1_E1E2P1P2/{file_name}.tsv",sep="\t")
    out_path = f"mutate_pairs_{file_name}.csv"
    write_pair_mutation(df_regions,seq_extractor,jaspar_annotator,device,out_path)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_name",type=str)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    analysis(args.file_name,args.device)
    logger.info(f"Done with {args.file_name}")