#!/usr/bin/python3

"""
This modules mutate bases in sequences and calculate ism_delta
Can be used in two mode: single and bulk
Single mode: Use in python script only: input one sequence, no intermediate directory is generated.
Bulk mode: Input sequences, create intermediate directory containing all mutated files for paralle processing
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




def saturation_SNPs_all_seqs(seqs,temp_out_path):
    """
    Write all mutated sequences of all sequences in seq_file into dir_temp
    Args:
        seq_file: a file containing sequences
        seq_colname: the column name of the sequences 
    Output:
        NULL
    """
    colnames=["mutated_sequence","mutated_location","Refseq_idx"]
    with open(temp_out_path) as f:
        f.write(",".join(colnames)+"\n")

    for i,seq in enumerate(seqs):
        df=pd.DataFrame({
                         "mutated_sequence":saturation_SNPs_one_seq(seq),
                         "mutated_location":[item for item in list(range(0,len(seq)))for _ in range(3)]
        })
        df["Refseq_idx"]="Ref"+str(i)
        df.to_csv(temp_out_path,index=False,header=False,mode="a")



def aggregate_mutation_prediction(key):
    """
    Calculate average effect of mutating each location
    """
    df_mutated_info=pd.read_csv("Temp_ism/sequences_mutated.csv")
    df_mutated_prediction=pd.read_csv("Temp_ism/prediction_mutated.csv")
    assert df_mutated_info.shape[0]==df_mutated_prediction.shape[0]
    df_mutated=pd.concat([df_mutated_info,df_mutated_prediction],axis=1)
    df_mutated=df_mutated.groupby(["Refseq_idx","mutated_location"])["signal_"+key].mean().reset_index()
    return df_mutated


def calculate_ism_delta(device,model,seq_file,seq_colname,key,out_path):
    """ Sequence needs to be 600 bp"""
    # make directory for temp result
    random_number = random.randint(100000, 999999)
    dir_temp=f"Temp_ism_{random_number}"
    
    if not os.path.isdir("Temp_ism/"):
        os.makedirs("Temp_ism/")

    # write signal for reference sequence
    write_predictions(device,model,seq_file,
                            seq_colname,"Temp_ism/prediction_reference.csv",key,
                            variable_length=True,batch_size=8192)

    # write all mutated sequences
    seqs=pd.read_csv(seq_file).loc[:,seq_colname]
    saturation_SNPs_all_seqs(seqs)

    # write signal for mutated sequences
    write_predictions(device,model,"Temp_ism/sequences_mutated.csv",
                            "mutated_sequence","Temp_ism/prediction_mutated.csv",key,
                            variable_length=True,batch_size=8192)

    # aggretate sequence signal
    prediction_mutated=aggregate_mutation_prediction(key)
    prediction_reference=pd.read_csv("Temp_ism/prediction_reference.csv")
    refseq_idx=prediction_mutated["Refseq_idx"].str.replace('Ref', '').astype(int).values
    prediction_mutated["reference_signal"]=prediction_reference["signal_"+key].values[refseq_idx]
    prediction_mutated["ism_delta"]=prediction_mutated["reference_signal"]-prediction_mutated["signal_"+key]
    ism=prediction_mutated.loc[:,"ism_delta"].values.reshape(-1,600)
    pd.DataFrame(ism).to_csv(out_path,index=False)

    shutil.rmtree("Temp_ism/")





