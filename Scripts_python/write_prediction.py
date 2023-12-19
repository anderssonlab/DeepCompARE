#!/usr/bin/python3

"""
This module is both: 
a commandline tool: Given sequence file, write model prediction to csv.
a python function: compute_predictions()
"""


from utils import encode, find_available_gpu
import json
import numpy as np
import pandas as pd
import torch
import os
import argparse

index_dict=json.load(open("/isdata/alab/people/pcr980/DeepCompare/Scripts_python/index_dict.json"))



def _crop_seq(seq,padding="both_ends"):
    if len(seq)==600:
        return seq
    if len(seq)<=600:
        if padding=="both_ends":
            left_length=np.floor((600-len(seq))/2).astype(int)
            right_length=np.ceil((600-len(seq))/2).astype(int)
            return "N"*left_length+seq+"N"*right_length
        elif padding=="from_left":
            return "N"*(600-len(seq))+seq
        elif padding=="from_right":
            return seq+"N"*(600-len(seq))
        else:
            raise ValueError("Parameter padding unrecognized.Please enter 'both_ends', 'from_left' or 'from_right'")
    # when seq is longer than 600"
    start_idx=(len(seq)-600)//2
    return seq[start_idx:start_idx+600]


def _seq2x(seqs,device):
    if isinstance(seqs,str):
        seqs=[seqs]
    X=np.zeros((len(seqs),len(seqs[0]),4))
    X=np.array(list(map(encode,seqs))).transpose(0, 2, 1)
    if not device:
        return X
    X=torch.tensor(X,device=device).float()
    return X







def compute_predictions(seqs,model=False,device=False):
    """
    Aimed to be used in other python scripts, as a part of bigger analyses
    Args:
        seqs: a list of strings, shouldn't contain more than 10k sequences
        model: optional, loaded torch model. If absent, just use default model.
        device: optional, a torch.device(). If absent, find first GPU that has more than 5GB memory.
    Output: 
        numpy array of model output
    """
    if not device:
        gpu_idx=find_available_gpu()
        device=torch.device(f"cuda:{gpu_idx}")
    if not model:
        model=torch.load("/isdata/alab/people/pcr980/DeepCompare/AstigCRConv5D_Dataset_final_CR_Trainer/model.h5",map_location=device)
        model.eval()
    X=_seq2x(seqs,device)
    with torch.no_grad():
        res=model(X).cpu().detach().numpy()
    return res



def write_predictions(data_path,seq_colname,out_path,
                      gpu_idx="infer",
                      model_dir="/isdata/alab/people/pcr980/DeepCompare/AstigCRConv5D_Dataset_final_CR_Trainer/",
                      variable_length=False,
                      batch_size=16384):
    """
    Aimed to be used as a independent commandline tool. 
    Args:
        gpu_idx: optional, gpu index, from "0" to "7", or a torch.device()
        model_dir: optional, need to be absolute path of the directory that stores "model.h5"
    Output:
        NULL, but write model output to csv file
    """
    model=torch.load(os.path.join(model_dir, "model.h5"))
    model.eval()
    if gpu_idx=="infer":
        gpu_idx=find_available_gpu()
    device=torch.device(f"cuda:{gpu_idx}")
    model.to(device)
    for chunk in pd.read_csv(data_path,chunksize=batch_size,header=0):
        query_seqs=list(chunk.loc[:,seq_colname])
        if variable_length:
            query_seqs=list(map(_crop_seq,query_seqs))
        df=pd.DataFrame(compute_predictions(query_seqs,model,device))
        df.to_csv(out_path, mode='a',index=False,header=False)
        



# export functions
__all__ = ['compute_predictions',
           'write_predictions']

if __name__=="__main__":
    # parse parameters
    parser = argparse.ArgumentParser(description="Input parameters for prediction")
    parser.add_argument("-d", "--data_path", type=str, required=True, help="Data path")
    parser.add_argument("-c", "--seq_colname", type=str, required=True, help="Name of column that contains the sequences")
    parser.add_argument("-o", "--out_path", type=str, help="Output path")
    parser.add_argument("-m", "--model", type=str, default="/isdata/alab/people/pcr980/DeepCompare/AstigCRConv5D_Dataset_final_CR_Trainer/",help="Model name")
    parser.add_argument("-g", "--gpu_index", type=str, default="infer",help="GPU index")
    parser.add_argument("-v", "--variable_length", type=bool, default=False, help="Is input sequence of variable lengths? (default False)")
    parser.add_argument("-b", "--batch_size", type=int, default=8192, help="Batch size (default 8192)")
    
    args=parser.parse_args()
    
    if args.gpu_index=="infer":
        args.gpu_idx=find_available_gpu()
    
    write_predictions(args.data_path, args.seq_colname, args.out_path, 
                      args.gpu_index, args.model, args.variable_length, args.batch_size)


