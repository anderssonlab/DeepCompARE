#!/usr/bin/python3
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import seq2x, find_available_gpu, SeqExtractor

import pandas as pd
import numpy as np
import torch
from captum.attr import InputXGradient
from kipoiseq import Interval
import re



def backprop(interpreter, x, target):
    return interpreter.attribute(x, target=target).cpu().detach().numpy().sum(axis=1).squeeze()



def _extract_numbers(s):
    """
    Solely for sorting row names
    """
    return list(map(int, re.findall(r'\d+', s)))



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
        model: if False, load default model
        device: if False, load gpu found by find_available_gpu()
        start_idx: if seqs is a list, start_idx is the index of the first sequence.
                   if function is called alone, start_idx is 0.
                   if function is called by write_gradxinp_from_seq(), start_idx is passed by write_gradxinp_from_seq().
    Return:
        attributions: attributions in pandas dataframe, with "SeqX_trackY" as row names.
                      SeqX has same order as seqs
    """
    # to do: enable bulk (batch) mode
    
    # add default model and device
    model.to(device)
    interpreter = InputXGradient(model)
    
    # convert seqs to x
    x=seq2x(seqs,device)
    x.requires_grad_()
    
    # if single target, convert to list
    if isinstance(targets,int):
        targets=[targets]

    attributions=[]
    indices=[]
    for target in targets:
        attributions.append(backprop(interpreter, x, target))
        indices+=[f"Seq{seq_idx}_Track{target}" for seq_idx in range(start_idx,start_idx+len(seqs))]
    
    attributions=np.array(attributions).squeeze()
    attributions=attributions.reshape((len(targets)*len(seqs),len(seqs[0])))
    df=pd.DataFrame(attributions,index=indices)
    
    # reorder rows by indices
    sorted_indices = sorted(indices, key=lambda x: _extract_numbers(x))
    df=df.reindex(sorted_indices)
    return df
    




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
    
    

    
    
    
def write_gradxinp_from_bed_file(bed_file,out_path,targets=list(range(16)),
                                 model=torch.load("/isdata/alab/people/pcr980/DeepCompare/DeepCompare_model/model.h5"),
                                 device=torch.device("cuda:"+find_available_gpu()),
                                 batch_size=4096):
    extractor = SeqExtractor()
    bed_df=pd.read_csv(bed_file,header=None,sep="\t")
    # resize to 600 bp
    
    intervals = [Interval(chrom,start,end) for chrom,start,end in zip(bed_df.iloc[:,0],bed_df.iloc[:,1],bed_df.iloc[:,2])]
    seqs = [extractor.extract(interval) for interval in intervals]
    write_gradxinp_from_seq(seqs,out_path,targets=targets,model=model,device=device,batch_size=batch_size)




def read_gradxinp(gradxinp_file,track_num):
    """
    Read gradxinp from gradxinp_file, subset by track_num
    """
    gradxinp_df=pd.read_csv(gradxinp_file,header=None,index_col=0)
    # Given that indices are composed of "SeqX_TrackY", we can subset to contain only "_Track{track_num}"
    gradxinp_df=gradxinp_df[gradxinp_df.index.str.contains(f"_Track{track_num}$")]
    return gradxinp_df