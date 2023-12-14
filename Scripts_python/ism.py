#!/usr/bin/python3

"""
This modules mutate bases in sequences and calculate ism_delta
Can be used in two mode: single and bulk
Single mode: Use in python script only: input one sequence, no intermediate directory is generated.
Bulk mode: Input sequences, create intermediate directory containing all mutated files for paralle processing
Files created in Bulk mode:
    sequences_mutated.csv
    prediction_mutated_seqs.csv
    prediction_ref_seqs.csv
"""


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from write_prediction import write_predictions

import random
import shutil
import os
import pandas as pd



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
        seq_name: optional
        out_path: option
    Output:
        mutated sequences in a list
    """
    mutated_seqs_nested=[_mutate_to_3_alt_bases(seq,location) for location in range(len(seq))]
    mutated_seqs_flat=[item for sublist in mutated_seqs_nested for item in sublist]
    return mutated_seqs_flat




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
    out_path=os.path.join(dir_temp,"sequences_mutated.csv")
    with open(out_path) as f:
        f.write(",".join(colnames)+"\n")

    for i,seq in enumerate(seqs):
        df=pd.DataFrame({
                         "mutated_sequence":saturation_SNPs_one_seq(seq),
                         "mutated_location":[item for item in list(range(0,len(seq)))for _ in range(3)]
        })
        df["Refseq_idx"]="Ref"+str(i)
        df.to_csv(out_path,index=False,header=False,mode="a")







def aggregate_mutation_prediction(dir_temp):
    """
    Calculate average effect of mutating each location
    """
    df_mutated_info=pd.read_csv("Temp_ism/sequences_mutated.csv")
    df_mutated_prediction=pd.read_csv("Temp_ism/prediction_mutated_seqs.csv")
    assert df_mutated_info.shape[0]==df_mutated_prediction.shape[0]
    df_mutated=pd.concat([df_mutated_info,df_mutated_prediction],axis=1)
    df_mutated=df_mutated.groupby(["Refseq_idx","mutated_location"])["signal_"+key].mean().reset_index()
    return df_mutated


def _create_dir_name():
    random_number = random.randint(100000, 999999)
    return f"Temp_ism_{random_number}"

def calculate_ism_delta(device,model_dir,seq_file,seq_colname,out_path):
    """ 
    Args:
        
    """
    # make directory for temp result
    dir_temp=_create_dir_name()
    
    while os.path.isdir(dir_temp):
        dir_temp=_create_dir_name()
    
    print(f"Create directory {dir_temp}")
    os.makedirs(dir_temp)

    # write signal for reference sequence
    write_predictions(device,model_dir,seq_file,seq_colname,
                      os.path.join(dir_temp,"prediction_ref_seqs.csv"))

    # write all mutated sequences
    seqs=pd.read_csv(seq_file).loc[:,seq_colname]
    write_saturation_SNPs_for_all_seqs(seqs,dir_temp)

    # write signal for mutated sequences
    write_predictions(device,model_dir,
                      os.path.join(dir_temp,"sequences_mutated.csv"),
                      "mutated_sequence",
                      os.path.join(dir_temp,"prediction_mutated_seqs.csv"))

    # aggretate sequence signal
    prediction_mutated=aggregate_mutation_prediction(dir_temp)
    prediction_reference=pd.read_csv(os.path.join(dir_temp,"prediction_ref_seqs.csv"))
    refseq_idx=prediction_mutated["Refseq_idx"].str.replace('Ref', '').astype(int).values
    prediction_mutated["reference_signal"]=prediction_reference["signal_"+key].values[refseq_idx]
    prediction_mutated["ism_delta"]=prediction_mutated["reference_signal"]-prediction_mutated["signal_"+key]
    ism=prediction_mutated.loc[:,"ism_delta"].values.reshape(-1,600)
    pd.DataFrame(ism).to_csv(out_path,index=False)

    shutil.rmtree(dir_temp)





