#!/usr/bin/python3

import torch
import argparse
import pandas as pd
from utils import find_available_gpu
from gradxinp import compute_gradxinp_from_seq




def write_gradxinp_from_seq(seqs,out_path,targets=list(range(16)),
                            model=torch.load("/isdata/alab/people/pcr980/DeepCompare/DeepCompare_model/model.h5"),
                            device=torch.device("cuda:"+find_available_gpu()),
                            batch_size=4096):
    """
    Aim to perform computation for large number of sequences, in batches
    Write attributions to out_path
    """
    model.to(device)
    start_idx=0
    for i in range(0,len(seqs),batch_size):
        chunk=seqs[i:i+batch_size]
        importance = compute_gradxinp_from_seq(chunk,targets=targets,model=model,device=device,start_idx=start_idx)
        importance.to_csv(out_path,mode="a",header=False)
        start_idx+=len(seqs)



def write_gradxinp_from_seq_file(seq_file,seq_colname,out_path,targets=list(range(16)),
                            model=torch.load("/isdata/alab/people/pcr980/DeepCompare/DeepCompare_model/model.h5"),
                            device=torch.device("cuda:"+find_available_gpu()),
                            batch_size=4096):
    model.to(device)
    
    start_idx=0
    for chunk in pd.read_csv(seq_file,chunksize=batch_size,header=0):
        seqs=chunk.loc[:,seq_colname].tolist()
        importance = compute_gradxinp_from_seq(seqs,targets=targets,model=model,device=device,start_idx=start_idx)
        importance.to_csv(out_path,mode="a",header=False)
        start_idx+=len(seqs)
        
        

if __name__=="__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--seq_file",type=str,help="path to sequence file")
    argparser.add_argument("--seq_colname",type=str,help="column name of sequence")
    argparser.add_argument("--out_path",type=str,help="path to output file")

    args=argparser.parse_args()
    write_gradxinp_from_seq_file(args.seq_file,args.seq_colname,args.out_path)
    