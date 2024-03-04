#!/usr/bin/python3

"""
This modules mutate bases in sequences and calculate ism_delta.
Sequences must have same length
The output format fits gradxinp.py
Can be used in two mode: single and bulk
"""
import random
import shutil
import os
import pandas as pd
import numpy as np
import argparse
import torch
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from utils import find_available_gpu, extract_numbers
from prediction import compute_predictions,write_predictions




#------------------------------
# Calculate ISM
#-----------------------------


def _mutate_single_base(seq,location):
    return seq[:location]+"N"+seq[location+1:]


def _mutate_one_seq(seq):
    mutated_seqs=[_mutate_single_base(seq,location) for location in range(len(seq))]
    return mutated_seqs


def mutate_all_seqs(seqs):
    mutated_seqs=[_mutate_one_seq(seq) for seq in seqs]
    mutated_seqs=[item for sublist in mutated_seqs for item in sublist]
    df_seqs_mut=pd.DataFrame({"mutated_sequence":mutated_seqs})
    df_seqs_mut.index=[f'Seq{i}_{j}' for i in range(len(seqs)) for j in range(len(seqs[0]))]
    return df_seqs_mut
    
    
def calculate_delta_score(pred_ref:pd.DataFrame,
                          pred_mut:pd.DataFrame,
                          seq_len=600)->pd.DataFrame:
    pred_ref = np.repeat(pred_ref.values, seq_len, axis=0)
    assert pred_mut.shape==pred_ref.shape
    deltas=pd.DataFrame(pred_ref-pred_mut) # so that positively contributing bases get positive values
    deltas.columns=[f"Track{i}" for i in range(16)]
    split_df = pred_mut.index.to_series().str.split('_', expand=True)
    deltas['Sequence'] = split_df[0]
    deltas['Location'] = split_df[1].astype(int)
    deltas.set_index(['Sequence', 'Location'], inplace=True)
    deltas = deltas.stack().unstack(level=1)
    deltas.index= [f'{seq}_{track}' for seq, track in deltas.index]
    # reorder rows by number
    sorted_indices = sorted(deltas.index, key=lambda x: extract_numbers(x))
    deltas=deltas.reindex(sorted_indices)
    return deltas


def compute_ism_score(seqs):
    """
    Calculate average effect of mutating each location
    """
    if isinstance(seqs,str):
        seqs=[seqs]
    df_seqs_mut=mutate_all_seqs(seqs)
    pred_ref = pd.DataFrame(compute_predictions(seqs))
    pred_mut = pd.DataFrame(compute_predictions(df_seqs_mut.mutated_sequence))
    pred_mut.index=df_seqs_mut.index
    df_ism=calculate_delta_score(pred_ref,pred_mut)
    return df_ism
    

def _create_dir_name():
    random_number = random.randint(1, 999999)
    return f"Temp_ism_{random_number}"

def create_dir():
    wd=os.getcwd()
    dir_temp=_create_dir_name()
    while os.path.isdir(os.path.join(wd,dir_temp)):
        dir_temp=_create_dir_name()
    logger.info(f"Create directory {dir_temp}")
    os.makedirs(os.path.join(wd,dir_temp))
    return os.path.join(wd,dir_temp)


def write_ism_score(seq_file,seq_colname,out_path, dir_temp=None,
                    device=torch.device("cuda:"+find_available_gpu())):

    # make directory for temp result, if dir_temp is not provided
    # if called from write_ism_N_bed, dir_temp is provided
    if dir_temp is None:
        dir_temp=create_dir()

    # write signal for reference sequence
    write_predictions(seq_file,seq_colname,os.path.join(dir_temp,"pred_ref.csv"),
                      device=device)

    # write all mutated sequences
    seqs=pd.read_csv(seq_file).loc[:,seq_colname]
    df_seqs_mut=mutate_all_seqs(seqs)
    df_seqs_mut.to_csv(os.path.join(dir_temp,"seqs_mut.csv"))

    # write signal for mutated sequences
    write_predictions(os.path.join(dir_temp,"seqs_mut.csv"),
        "mutated_sequence",
        os.path.join(dir_temp,"pred_mut.csv"),
        device=device)

    pred_ref=pd.read_csv(os.path.join(dir_temp,"pred_ref.csv"),header=None)
    pred_mut=pd.read_csv(os.path.join(dir_temp,"pred_mut.csv"),header=None)
    pred_mut.index=df_seqs_mut.index
    df_ism=calculate_delta_score(pred_ref,pred_mut)
    pd.DataFrame(df_ism).to_csv(out_path,header=False)
    shutil.rmtree(dir_temp)
    logger.info(f"Done with {seq_file}. Delete directory {dir_temp}")
    
    
        
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Calculate ISM')
    parser.add_argument('--seq_file', type=str, help='path to sequence file')
    parser.add_argument('--seq_colname', type=str, help='column name of sequences')
    parser.add_argument('--out_path', type=str, help='path to output file')
    parser.add_argument('--gpu', type=str, help='gpu id')
    args = parser.parse_args()
    
    write_ism_score(args.seq_file,args.seq_colname,args.out_path,
                    device=torch.device(f"cuda:{args.gpu}"))