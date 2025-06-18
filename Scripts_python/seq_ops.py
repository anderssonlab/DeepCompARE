
import numpy as np
import pandas as pd
import torch
from pyfaidx import Fasta
from collections import Counter
from loguru import logger

class SeqExtractor:
    """
    Convert bed regions to fasta sequences
    """
    def __init__(self,fasta_path):
        self.fasta = Fasta(fasta_path)
    #
    def get_seq(self,region):
        # deal with iterables
        if isinstance(region, list) or isinstance(region, tuple) and len(region)>2:
            chr=region[0]
            start=int(region[1])
            end=int(region[2])
        elif isinstance(region, str):
            # deal with format like chr:start-end
            chr, start_end = region.split(":")
            start, end = start_end.split("-")
            start, end = int(start), int(end)
        else:
            raise ValueError("Invalid input.")
        try:
            return str(self.fasta.get_seq(chr,start,end).seq).upper()
        except:
            logger.error(f"Region {chr} {start} {end} not found in the fasta file.")
            return None
            





def generate_random_seq(length):
    return ''.join(np.random.choice(["A", "C", "G", "T", "N"], length))  



    
def generate_random_seqs(num_seqs, length):
    return [generate_random_seq(length) for _ in range(num_seqs)]    
    



def ablate_motifs(seq, motif_starts, motif_ends):
    """
    Scramble the sequence between multiple motif start and end positions.
    Args:
        seq: A string of sequence.
        motif_starts: A list of integers for motif starts.
        motif_ends: A list of integers for motif ends.
    Returns:
        A string of scrambled sequence.
    """
    if isinstance(motif_starts, int):
        motif_starts = [motif_starts]
    if isinstance(motif_ends, int):
        motif_ends = [motif_ends]
    if len(motif_starts) != len(motif_ends):
        raise ValueError("motif_starts and motif_ends must have the same length")
    # Sort the motifs by start position
    motifs = sorted(zip(motif_starts, motif_ends), key=lambda x: x[0])
    # Initialize variables
    seq_ablated = ''
    previous_end = 0
    # Iterate and ablate each motif
    for start, end in motifs:
        if start < previous_end:
            raise ValueError("Overlapping motifs detected")
        end = end + 1  
        motif = seq[start:end]
        motif_scrambled = "N" * len(motif)  
        # Append non-motif and scrambled motif parts
        seq_ablated += seq[previous_end:start] + motif_scrambled
        previous_end = end
    # Append the remaining part of the sequence if any
    seq_ablated += seq[previous_end:]
    return seq_ablated



def encode(my_seq):
    mapping= { "A": [1, 0, 0, 0], "C": [0, 1, 0, 0],"G": [0, 0, 1, 0],"T": [0, 0, 0, 1],"N": [0, 0, 0, 0]}
    return np.array([mapping[i] for i in my_seq])




def dinucleotide_frequency(seq):
    """
    Args:
        seq: a string of DNA sequence
    Returns:
        a dictionary of dinucleotide frequency
    """
    dinucleotides = [seq[i:i+2] for i in range(len(seq)-1)]
    dinucleotide_freq=Counter(dinucleotides)
    all_dinucleotides = [a+b for a in 'ACGT' for b in 'ACGT']
    complete_freq = {dn: dinucleotide_freq.get(dn, 0) for dn in all_dinucleotides}
    df = pd.DataFrame([complete_freq])
    return df




def dinucleotide_frequencies(seqs):
    """
    Args:
        seqs: a list of strings, or a single string
    Returns:
        a dataframe of dinucleotide frequencies
    """
    if isinstance(seqs,str):
        seqs=[seqs]
    if hasattr(seqs,"values"):
        seqs=seqs.values
    dinucleotide_freqs = [dinucleotide_frequency(seq) for seq in seqs]
    return pd.concat(dinucleotide_freqs, ignore_index=True)




def trinucleotide_frequency(seq):
    trinucleotides = [seq[i:i+3] for i in range(len(seq)-2)]
    trinucleotide_freq=Counter(trinucleotides)
    all_trinucleotides = [a+b+c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']
    complete_freq = {tn: trinucleotide_freq.get(tn, 0) for tn in all_trinucleotides}
    df = pd.DataFrame([complete_freq])
    return df




def trinucleotide_frequencies(seqs):
    """
    Args:
        seqs: a list of strings, or a single string
    Returns:
        a dataframe of trinucleotide frequencies
    """
    if isinstance(seqs,str):
        seqs=[seqs]
    if hasattr(seqs,"values"):
        seqs=seqs.values
    trinucleotide_freqs = [trinucleotide_frequency(seq) for seq in seqs]
    return pd.concat(trinucleotide_freqs, ignore_index=True)
    



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


def resize_seq(seq,padding="both_ends"):
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





def shuffle_dinucleotides(seq):
    evens = seq[0::2]
    odds = seq[1::2]
    zipped_pairs = zip(odds, evens)
    shuffled_sequence = ''.join(a + b for a, b in zipped_pairs)
    # Handle the case where the original sequence length is odd
    if len(seq) % 2 != 0:
        shuffled_sequence += seq[-1]
    return shuffled_sequence

