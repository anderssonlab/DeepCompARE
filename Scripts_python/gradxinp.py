"""
This modules calculates gradxinp.
Output format: 
    Row names are "SeqX_TrackY", where X is the index of sequence, Y is the index of track.
    

"""


from utils import find_available_gpu,extract_numbers
from seq_ops import seq2x

import pandas as pd
import numpy as np
import torch
from captum.attr import InputXGradient




def backprop(interpreter, x, target):
    return interpreter.attribute(x, target=target).cpu().detach().numpy().sum(axis=1).squeeze()




def compute_gradxinp_from_seq(seqs,targets=list(range(16)),
                              model=torch.load("/isdata/alab/people/pcr980/DeepCompare/DeepCompare_model/model.h5"),
                              device=torch.device("cuda:"+find_available_gpu()),
                              start_idx=0):
    """
    Aim to perform computation for small number of sequences
    Return attributions in pandas dataframe
    Args:
        seqs: sequence, either a string or a list of strings. Doesn't support batch operation, sow len(seq) should be small
        target: number of track to use. If False, return attributions for all targets
        start_idx: if seqs is a list, start_idx is the index of the first sequence.
                   if function is called alone, start_idx is 0.
                   if function is called by write_gradxinp_from_seq(), start_idx is passed by write_gradxinp_from_seq().
    Return:
        attributions: attributions in pandas dataframe, with "SeqX_trackY" as row names.
                      SeqX has same order as seqs
    """
    if isinstance(seqs,str):
        seqs=[seqs]
    if isinstance(targets,int):
        targets=[targets]
    
    model.to(device)
    interpreter = InputXGradient(model)
    x=seq2x(seqs,device)
    x.requires_grad_()
    attributions=[]
    indices=[]
    # to do: allow batch operation
    for target in targets:
        attributions.append(backprop(interpreter, x, target))
        indices+=[f"Seq{seq_idx}_Track{target}" for seq_idx in range(start_idx,start_idx+len(seqs))]
    attributions=np.array(attributions).squeeze()
    attributions=attributions.reshape((len(targets)*len(seqs),len(seqs[0])))
    df=pd.DataFrame(attributions,index=indices)
    # reorder rows by indices
    sorted_indices = sorted(indices, key=lambda x: extract_numbers(x))
    df=df.reindex(sorted_indices)
    return df
    

