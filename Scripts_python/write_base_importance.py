import torch
import argparse
import pandas as pd
from utils import find_available_gpu
from seq_ops import SeqExtractor
from gradxinp import  compute_gradxinp
from in_silico_mutagenesis import compute_mutagenesis_score





def write_base_imp(file_name,imp_type,out_path,hg="hg38",
                   device=torch.device("cuda:"+find_available_gpu())):
    if hg=="hg38":
        seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    elif hg=="hg19":
        seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg19.fa")
    else:
        raise ValueError("hg must be either hg38 or hg19")

    if imp_type=="gradxinp":
        df_imp=compute_gradxinp(file_name,
                                seq_extractor=seq_extractor,
                                device=device)
    elif imp_type=="ism":
        df_imp=compute_mutagenesis_score(file_name,
                                         seq_extractor=seq_extractor,
                                         mode="ism",
                                         device=device)
    elif imp_type=="isa":
        df_imp=compute_mutagenesis_score(file_name,
                                         seq_extractor=seq_extractor,
                                         mode="isa",
                                         device=device)
    else:
        raise ValueError("imp_type must be either gradxinp, ism or isa")
    df_imp.to_csv(out_path,header=False)






if __name__=="__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument("--file_name",type=str,help="path to sequence file")
    parser.add_argument("--imp_type",type=str,help="type of importance score")
    parser.add_argument("--out_path",type=str,help="path to output file")
    parser.add_argument("--device",type=str,help="gpu index")

    args=parser.parse_args()
    write_base_imp(args.file_name,
                   args.imp_type,
                   args.out_path,
                   device=torch.device(f"cuda:{args.device}"))
