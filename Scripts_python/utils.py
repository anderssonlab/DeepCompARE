import numpy as np
import pandas as pd
import torch
import pynvml
import re

from pyfaidx import Fasta
from kipoiseq import Interval
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_random_seq(length):
    return ''.join(np.random.choice(["A", "C", "G", "T", "N"], length))  
    
def generate_random_seqs(num_seqs, length):
    return [generate_random_seq(length) for _ in range(num_seqs)]    
    

def encode(my_seq):
    mapping= { "A": [1, 0, 0, 0], "C": [0, 1, 0, 0],"G": [0, 0, 1, 0],"T": [0, 0, 0, 1],"N": [0, 0, 0, 0]}
    return np.array([mapping[i] for i in my_seq])


def shift_seq(seq,direction,dist):
    """ 
    add "N"*dist to object seq from directin 
    """
    if direction=="right": # add 0 from left,subset from left
        return ("N"*dist+seq)[0:600]
    elif direction=="left": # add 0 from right, subset from right
        return (seq+"N"*dist)[dist:dist+600]
    else:
        raise ValueError("parameter direction not recognized!")




def seq2x(seqs,device=False):
    """
    Args:
        seqs: a list of strings, or a single string
        device: a torch.device()
    Returns:
        numpy array of shape (len(seqs),4,len(seqs[0]))
    """
    if isinstance(seqs,str):
        seqs=[seqs]
    if hasattr(seqs,"values"):
        seqs=seqs.values
    X=np.zeros((len(seqs),len(seqs[0]),4))
    X=np.array(list(map(encode,seqs))).transpose(0, 2, 1)
    if not device:
        return X
    X=torch.tensor(X,device=device).float()
    return X



def find_available_gpu():
    """
    Find first available GPU
    Returns:
        str: GPU id
    """
    pynvml.nvmlInit()
    deviceCount = pynvml.nvmlDeviceGetCount()
    for i in range(deviceCount):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        mem = pynvml.nvmlDeviceGetMemoryInfo(handle)
        if mem.free/1024**3>20:
            return str(i)
    raise ValueError("No available GPU found!")



class SeqExtractor:
    """
    Convert bed regions to fasta sequences
    copied from enformer_usage.ipynb
    """
    def __init__(self):
        self.fasta = Fasta("/binf-isilon/alab/people/pcr980/Resource/hg38.fa")
    
    def get_seq(self,chr,start,end):
        return str(self.fasta.get_seq(chr,start,end).seq).upper()


def extract_numbers(s):
    """
    Solely for sorting row names
    """
    return list(map(int, re.findall(r'\d+', s)))
    
    
def read_featimp(featimp_file,track_num):
    """
    Read featimp from featimp_file, subset by track_num
    Featimp is either gradxinp or ism, no header
    """
    featimp_df=pd.read_csv(featimp_file,header=None,index_col=0)
    # Given that indices are composed of "SeqX_TrackY", we can subset to contain only "_Track{track_num}"
    featimp_df=featimp_df[featimp_df.index.str.contains(f"_Track{track_num}$")]
    return featimp_df

    
def write_fasta(out_fname,importance,seqs_ids,seqs): # one number per location
    with open(out_fname,"a+") as f:
        for i in range(len(seqs)):
            seq_rec=SeqRecord(Seq(seqs[i]))
            seq_rec.id=seqs_ids[i]
            SeqIO.write(seq_rec,f,"fasta")
            f.write(' '.join(map(lambda x: f'{float(x):.5f}', importance[i])))
            f.write("\n")
