#!/usr/bin/python3

"""
This module is both: 
a commandline tool: Given sequence file, write model prediction to csv.
a python function: compute_predictions()
"""


from encode_sequence import encode
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


def _seq2x(seqs,device=False):
    X=np.zeros((len(seqs),len(seqs[0]),4))
    X=np.array(list(map(encode,seqs))).transpose(0, 2, 1)
    if not device:
        return X
    X=torch.tensor(X,device=device).float()
    return X








def compute_predictions(model,seqs,device):
    """
    Aimed to be used in other python script, as a part of bigger analyses
    Args:
        model: loaded torch model
        seqs: a list of strings
        device: a torch.device()
    Output: 
        numpy array of model output
    """
    X=_seq2x(seqs,device)
    with torch.no_grad():
        res=model(X).cpu().detach().numpy()
    return res


def write_predictions(gpu_idx,model_dir,data_path,seq_colname,out_path,variable_length=False,batch_size=8192):
    """
    Aimed to be used as a independent commandline tool. 
    Args:
        gpu_idx: gpu index, from "0" to "7"
        model_dir: need to be absolute path of the directory that stores "model.h5"
    Output:
        NULL, but write model output to csv file
    """
    model=torch.load(os.path.join(model_dir, "model.h5"))
    model.eval()
    device=torch.device(f"cuda:{gpu_idx}")
    model.to(device)
    for chunk in pd.read_csv(data_path,chunksize=batch_size,header=0):
        query_seqs=list(chunk.loc[:,seq_colname])
        if variable_length:
            query_seqs=list(map(_crop_seq,query_seqs))
        df=pd.DataFrame(compute_predictions(model,query_seqs,device))
        df.to_csv(out_path, mode='a',index=False,header=False)
        



# export functions
__all__ = ['compute_predictions',
           'write_predictions']

if __name__=="__main__":
    # parse parameters
    parser = argparse.ArgumentParser(description="Input parameters for prediction")
    parser.add_argument("-g", "--gpu_index", type=str, required=True, help="GPU index")
    parser.add_argument("-m", "--model", type=str, required=True, help="Model name")
    parser.add_argument("-d", "--data_path", type=str, required=True, help="Data path")
    parser.add_argument("-c", "--seq_colname", type=str, required=True, help="Name of column that contains the sequences")
    parser.add_argument("-o", "--out_path", type=str, help="Output path")
    parser.add_argument("-v", "--variable_length", type=bool, default=False, help="Is input sequence of variable lengths? (default False)")
    parser.add_argument("-b", "--batch_size", type=int, default=8192, help="Batch size (default 8192)")
    
    args=parser.parse_args()
    write_predictions(args.gpu_index, args.model, args.data_path, args.seq_colname, 
                      args.out_path, args.variable_length, args.batch_size)


