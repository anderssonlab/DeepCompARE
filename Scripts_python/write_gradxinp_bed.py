#!/usr/bin/python3

import torch
import argparse
import pandas as pd
from utils import find_available_gpu
from seq_ops import SeqExtractor
from write_gradxinp_seq import  write_gradxinp_from_seq



    
def write_gradxinp_from_bed_file(bed_file,out_path,targets=list(range(16)),
                                 model=torch.load("/isdata/alab/people/pcr980/DeepCompare/DeepCompare_model/model.h5"),
                                 device=torch.device("cuda:"+find_available_gpu()),
                                 batch_size=4096):
    extractor = SeqExtractor()
    bed_df=pd.read_csv(bed_file,header=None,sep="\t")
    intervals = [(chrom,start,end) for chrom,start,end in zip(bed_df.iloc[:,0],bed_df.iloc[:,1],bed_df.iloc[:,2]-1)] # -1 to not include the end
    seqs = [extractor.get_seq(*interval) for interval in intervals]
    write_gradxinp_from_seq(seqs,out_path,targets=targets,model=model,device=device,batch_size=batch_size)






if __name__=="__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--bed_file",type=str,help="path to sequence file")
    argparser.add_argument("--out_path",type=str,help="path to output file")
    argparser.add_argument("--device",type=str,help="gpu index")

    args=argparser.parse_args()
    write_gradxinp_from_bed_file(args.bed_file,
                                 args.out_path,
                                 device=torch.device(f"cuda:{args.device}"))
    