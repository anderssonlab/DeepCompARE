#!/usr/bin/python3

############################################################
import sys
sys.path.insert(1, "/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")


import os
import torch
import random
import shutil
import argparse
import random

import numpy as np
import pandas as pd

from typing import List
from loguru import logger

from utils import find_available_gpu
from seq_ops import resize_seq
from prediction import compute_predictions, write_predictions
from write_ism_N_seq import create_dir, calculate_delta_score

###############################################################################

BATCH_SIZE=8192


def _mutate_single_base(seq: str, location: int) -> List[str]:
    all_bases = ['A', 'T', 'C', 'G', 'N']
    
    all_bases.remove(seq[location])

    return [seq[:location] + base + seq[location+1:] for base in all_bases]


def _mutate_one_seq(seq: str) -> List[List[str]]:
    mutate_seqs = [[], [], [], []]
    for location in range(len(seq)):
        s1, s2, s3, s4 = _mutate_single_base(seq, location)
        mutate_seqs[0].append(s1)
        mutate_seqs[1].append(s2)
        mutate_seqs[2].append(s3)
        mutate_seqs[3].append(s4)
    return np.array(mutate_seqs)


def mutate_all_seqs(seqs: str) -> pd.DataFrame:
    all_mutate_seqs = []
    for seq in seqs:
        all_mutate_seqs.append(_mutate_one_seq(seq))
    
    all_mutate_seqs = np.hstack(all_mutate_seqs).T
    a1, a2 = np.meshgrid(np.arange(len(seqs)).astype('str'),
                            np.arange(len(seqs[0])).astype('str'))
    index = np.char.add(np.char.add(np.char.add("Seq", a1.T.ravel()), '_'), 
                        a2.T.ravel())
    
    return pd.DataFrame(data=all_mutate_seqs, index=index, 
                        columns=[f'mutated_sequence{i}' for i in range(1,5)])


def compute_ism_score(seqs):
    """
    Calculate average effect of mutating each location
    """
    if isinstance(seqs,str):
        seqs=[seqs]
    df_seqs_mut = mutate_all_seqs(seqs)
    pred_ref = pd.DataFrame(compute_predictions(seqs))

    pred_mut = compute_predictions(df_seqs_mut.mutated_sequence1)
    pred_mut += compute_predictions(df_seqs_mut.mutated_sequence2)
    pred_mut += compute_predictions(df_seqs_mut.mutated_sequence3)
    pred_mut += compute_predictions(df_seqs_mut.mutated_sequence4)
    pred_mut = pd.DataFrame(pred_mut/4)
    pred_mut.index=df_seqs_mut.index
    
    df_ism = calculate_delta_score(pred_ref,pred_mut)
    return df_ism


def write_A_predictions(data_path: os.PathLike,
                      seq_colnames: list,
                      out_path: os.PathLike,
                      model=torch.load("/isdata/alab/people/pcr980/DeepCompare/DeepCompare_model/model.h5"),
                      device=torch.device("cuda:"+find_available_gpu()),
                      variable_length=False,
                      batch_size=BATCH_SIZE):
    """
    Aimed to be used as a independent commandline tool. 
    Args:
        gpu_idx: optional, gpu index, from "0" to "7", or a torch.device()
        model_dir: optional, need to be absolute path of the directory that stores "model.h5"
    Output:
        NULL, but write model output to csv file. No header
    """
    model.eval()
    model.to(device)
    for chunk in pd.read_csv(data_path,chunksize=batch_size,header=0):
        a = []
        for seq_colname in seq_colnames:
            query_seqs=list(chunk.loc[:,seq_colname])
            if variable_length:
                query_seqs=list(map(resize_seq,query_seqs))
            
            computed_preds = compute_predictions(query_seqs, model, device)
            if len(a) == 0:
                a.append(computed_preds)
            else:
                a[0] += computed_preds
        df = pd.DataFrame(a[0]/4)
        df.to_csv(out_path, mode='a',index=False,header=False)



def write_ism_score(seq_file,seq_colname,out_path, dir_temp=None,
                    device=torch.device("cuda:"+ find_available_gpu())):
 
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
    df_seqs_mut.to_csv(os.path.join(dir_temp, "seqs_mut.csv"))

    # write signal for mutated sequences
    write_A_predictions(os.path.join(dir_temp,"seqs_mut.csv"),
        [f"mutated_sequence{i}" for i in range(1,5)],
        os.path.join(dir_temp,"pred_mut.csv"),
        device=device)

    pred_ref=pd.read_csv(os.path.join(dir_temp,"pred_ref.csv"), header=None)
    pred_mut=pd.read_csv(os.path.join(dir_temp,"pred_mut.csv"), header=None)
    pred_mut.index=df_seqs_mut.index
    df_ism=calculate_delta_score(pred_ref,pred_mut)
    pd.DataFrame(df_ism).to_csv(out_path,header=False)
    shutil.rmtree(dir_temp)
    logger.info(f"Done with {seq_file}. Delete directory {dir_temp}")