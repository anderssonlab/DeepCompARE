#!/usr/bin/python3

"""
This modules mutate bases in sequences and calculate ism_delta
Can be used in two mode: single and bulk
Single mode: Use in python script only: input one sequence, no intermediate directory is generated.
Bulk mode: Input sequences, create intermediate directory containing all mutated files for parallel processing
Files created in Bulk mode:
    seqs_mutated.csv
    pred_mutated.csv
    pred_ref.csv
"""
import random
import shutil
import os
import pandas as pd

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from prediction import compute_predictions,write_predictions




#------------------------------
# Calculate ISM
#------------------------------

def _mutate_single_base(seq,alt_base,location):
    """
    Args:
        seq: a string with capital letter
        alt_base: a character
        location: index of the base to change to alt_base, location is zero-based
    Output:
        a sequence of type string, with base at location chaged to alt_base
    """
    base_list = list(seq)
    base_list[location]=alt_base
    return ''.join(base_list)



def _mutate_to_3_alt_bases(seq,location):
    """
    Args:
        seq: a string with capital letter
        location: index of the base to change to alt_base, location is zero-based
    Output: 
        a list containing 3 sequences. At location, each sequence is mutated to one alternative bases
    """
    alt_bases=["A","C","G","T"]
    alt_bases.remove(seq[location])
    return [_mutate_single_base(seq,alt_base,location) for alt_base in alt_bases]



def saturation_SNPs_one_seq(seq):
    """
    Args:
        seq: one sequence
    Output:
        mutated sequences in a list
    """
    mutated_seqs_nested=[_mutate_to_3_alt_bases(seq,location) for location in range(len(seq))]
    mutated_seqs_flat=[item for sublist in mutated_seqs_nested for item in sublist]
    df=pd.DataFrame({"mutated_sequence":mutated_seqs_flat,
                     "mutated_location":[item for item in list(range(0,len(seq)))for _ in range(3)]})
    return df




def write_saturation_SNPs_for_all_seqs(seqs,dir_temp):
    """
    Write all mutated sequences of all sequences in seq_file into dir_temp
    Args:
        seq_file: a file containing sequences
        seq_colname: the column name of the sequences 
    Output:
        NULL, but write mutated sequences to csv file.
    """
    
    colnames=["mutated_sequence","mutated_location","Refseq_idx"]
    out_path=os.path.join(dir_temp,"seqs_mutated.csv")
    with open(out_path,"w") as f:
        f.write(",".join(colnames)+"\n")

    for i,seq in enumerate(seqs):
        df=saturation_SNPs_one_seq(seq)
        df["Refseq_idx"]="Ref"+str(i)
        df.to_csv(out_path,index=False,header=False,mode="a")





def _aggregate_mutation_prediction(seq=False,dir_temp=False):
    """
    Calculate average effect of mutating each location
    """
    if seq:
        assert all(char in ['A', 'C', 'G', 'T', 'N'] for char in seq)
        df_mutated_seqs=saturation_SNPs_one_seq(seq)
        pred_ref = compute_predictions(seq)
        pred_mutated = pd.DataFrame(compute_predictions(df_mutated_seqs.mutated_sequence))
        deltas=pred_mutated-pred_ref # broadcasting
        df_mutated=pd.concat([df_mutated_seqs,deltas],axis=1)
        df_ism=df_mutated.groupby(["mutated_location"]).mean()
    
    if dir_temp:
        df_mutated_seqs=pd.read_csv(os.path.join(dir_temp,"seqs_mutated.csv"))
        pred_mutated=pd.read_csv(os.path.join(dir_temp,"pred_mutated.csv"),header=None)
        assert df_mutated_seqs.shape[0]==pred_mutated.shape[0]
        pred_ref=pd.read_csv(os.path.join(dir_temp,"pred_ref.csv"),header=None)
        pred_ref["Refseq_idx"]= ["Ref" + str(i) for i in range(pred_ref.shape[0])]
        pred_ref.set_index("Refseq_idx",inplace=True)
        pred_mutated = pred_mutated.subtract(pred_ref.loc[df_mutated_seqs["Refseq_idx"]].values)
        df_ism=pd.concat([df_mutated_seqs,pred_mutated],axis=1)
        df_ism=df_ism.groupby(["Refseq_idx","mutated_location"]).mean().reset_index()

    return df_ism
        

def _create_dir_name():
    random_number = random.randint(1, 999999)
    return f"Temp_ism_{random_number}"


def calculate_ism_delta(*args):
    """ 
    """
    if len([*args])==3:
        seq_file,seq_colname,out_path=args
    
        # make directory for temp result
        dir_temp=_create_dir_name()
    
        while os.path.isdir(dir_temp):
            dir_temp=_create_dir_name()
    
        print(f"Create directory {dir_temp}")
        os.makedirs(dir_temp)

        # write signal for reference sequence
        write_predictions(seq_file,seq_colname,os.path.join(dir_temp,"pred_ref.csv"))

        # write all mutated sequences
        seqs=pd.read_csv(seq_file).loc[:,seq_colname]
        write_saturation_SNPs_for_all_seqs(seqs,dir_temp)

        # write signal for mutated sequences
        write_predictions(os.path.join(dir_temp,"seqs_mutated.csv"),
                         "mutated_sequence",
                         os.path.join(dir_temp,"pred_mutated.csv"))

       # aggretate sequence signal
        df_ism=_aggregate_mutation_prediction(dir_temp=dir_temp)
        pd.DataFrame(df_ism).to_csv(out_path,index=False)
        shutil.rmtree(dir_temp)
        print(f"ISM calculation complete!\nRemove directory {dir_temp}")
        
    elif len([*args])==1:
        if isinstance(args[0],str):
            print("Calculate ISM for one sequence")
            return _aggregate_mutation_prediction(seq=args[0])
        elif isinstance(args[0],list):
            # to do
            pass
        else:
            raise ValueError("Wrong type of argument")
        
    else:
        raise ValueError("Wrong number of arguments")
        


