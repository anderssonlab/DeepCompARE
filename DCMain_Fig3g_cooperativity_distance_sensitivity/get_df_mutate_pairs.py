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

def analysis(file_name,device,sep="\t"):
    # load data and tools
    df_regions = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_name}.bed",sep=sep)
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    if "hepg2" in file_name:
        jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                        score_thresh=500,
                                        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv",
                                        #chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed"
                                        )
    elif "k562" in file_name:
        jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                        score_thresh=500,
                                        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv",
                                        #chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed"
                                        )
    else:
        raise ValueError("file_name should contain hepg2 or k562")
    #
    out_path = f"mutate_pairs_cnn6_{file_name}.csv"
    # TODO: change model for prediction
    write_pair_mutation(df_regions,seq_extractor,jaspar_annotator,device,out_path)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_name",type=str)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    analysis(args.file_name,args.device)
    logger.info(f"Done with {args.file_name}")