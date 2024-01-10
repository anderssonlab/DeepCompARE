import numpy as np
import torch

from pyfaidx import Fasta
from kipoiseq import Interval
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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




def seq2x(seqs,device):
    """
    Args:
        seqs: a list of strings, or a single string
        device: a torch.device()
    Returns:
        numpy array of shape (len(seqs),4,len(seqs[0]))
    """
    if isinstance(seqs,str):
        seqs=[seqs]
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
    if torch.cuda.is_available():
        # Number of GPUs available
        num_gpus = torch.cuda.device_count()
        for i in range(num_gpus):
            device = torch.device("cuda:"+str(i))
            total_memory = torch.cuda.get_device_properties(i).total_memory
            allocated_memory = torch.cuda.memory_reserved(device)
            free_memory = total_memory - allocated_memory
            
            if free_memory>(1024**3)*20:
                return str(i)
            
        raise Exception("No GPU available.")
    
    
    
    
    
class SeqExtractor:
    """
    Convert bed regions to fasta sequences
    copied from enformer_usage.ipynb
    """
    def __init__(self):
        self.fasta = Fasta("/binf-isilon/alab/people/pcr980/Resource/hg38.fa")
    
    def get_seq(self,chr,start,end):
        return str(self.fasta.get_seq(chr,start,end).seq).upper()

    
    
    
def write_fasta(out_fname,importance,seqs_ids,seqs): # one number per location
    with open(out_fname,"a+") as f:
        for i in range(len(seqs)):
            seq_rec=SeqRecord(Seq(seqs[i]))
            seq_rec.id=seqs_ids[i]
            SeqIO.write(seq_rec,f,"fasta")
            f.write(' '.join(map(lambda x: f'{float(x):.5f}', importance[i])))
            f.write("\n")

