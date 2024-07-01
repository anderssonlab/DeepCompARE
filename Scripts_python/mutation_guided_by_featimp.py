import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Utils")
from encode_sequence import encode
from Scripts_python.prediction import predict_from_seq,seq2x, write_predictions
from write_feat_attr import compute_AttrxInp,write_attr

import shutil
import torch
import os
import pandas as pd
import numpy as np



#-------------------------------------------------------------------
# Mutate multiple locations, starting from most important sequences 
#--------------------------------------------------------------------
def find_first_two_consecutive_negs(arr):
    neg_indices=np.where(arr<0)[0]
    if len(neg_indices)==0:
        return len(arr)
    for i in range(len(neg_indices)-1):
        if neg_indices[i+1]==neg_indices[i]+1:
            return neg_indices[i]
    return len(arr)


def mutate_multiple_locations(seq,alt_base,locations):
    """
    change LOCATION of seq to alt_base
    locations must be np.array,containing multiple locations
    alt_base can be either character or np array
    """
    base_array = np.array(list(seq))
    base_array[locations]=alt_base
    base_list=base_array.tolist()
    return ''.join(base_list)

def calculate_attrxinp(seq,attr_4_rows):
    seq_one_hot=seq2x([seq])[0]
    return (attr_4_rows.values*seq_one_hot).sum(axis=0)


def find_most_deleterious_bases(reshaped_attr):
    # get columnwise rank of each attribution in increasing order
    rank_df=reshaped_attr.rank(axis=0) 
    # get columnwise idx of lowest value
    most_deleterious_bases=rank_df.apply(lambda x: np.where(x == 1)[0], axis=0).iloc[0,:].tolist()
    # convert 0123 to bases
    conversion_dict={0:"A",1:"C",2:"G",3:"T"}
    return np.array([conversion_dict[num] for num in most_deleterious_bases])


def reshape_attr(row):
    """
    Convert attribution of shape 1x2400 to 4x600
    """
    reshaped_df = pd.DataFrame()

    # Extract values for each condition and create new rows
    values_mod_0 = row.iloc[::4].values
    values_mod_1 = row.iloc[1::4].values
    values_mod_2 = row.iloc[2::4].values
    values_mod_3 = row.iloc[3::4].values
    # Create new rows and append them to the reshaped DataFrame
    new_row_0 = pd.DataFrame([values_mod_0])
    new_row_1 = pd.DataFrame([values_mod_1])
    new_row_2 = pd.DataFrame([values_mod_2])
    new_row_3 = pd.DataFrame([values_mod_3])
    reshaped_df = reshaped_df.append([new_row_0, new_row_1, new_row_2, new_row_3], ignore_index=True)
    return reshaped_df


def incremental_mutated_sequences_for_one_seq(seq,feat_attr,method):
    """
    Get all incrementally mutated sequences for a sequence seq
    if method=="zero": feat_attr has shape (1,600)
    if method=="most_deliterious": feat_attr has shape (1,2400)
    """
    if method=="zero":
        attrxinp=feat_attr
    elif method=="most_deleterious":
        # reshape attribution to 4 rows
        attr_4_rows=reshape_attr(feat_attr)
        attrxinp=calculate_attrxinp(seq,attr_4_rows)
    else:
        raise ValueError("Method not recognized!")
    location_order=np.flip(np.argsort(attrxinp))
    location_order_list = [location_order[:i] for i in range(1, len(location_order) + 1)]
    if method=="zero":
        return [mutate_multiple_locations(seq,"N",location) for location in location_order_list]
    else:
        # find worst base at every location
        bases=find_most_deleterious_bases(attr_4_rows)
        return [mutate_multiple_locations(seq,bases[location],location) for location in location_order_list]



def write_incremental_mutated_sequences_for_one_seq(seq,feat_attr,seq_name,method):
    """
    Write all incrementally mutated sequences for seq to csv
    """
    df=pd.DataFrame({
        "mutated_sequence":incremental_mutated_sequences_for_one_seq(seq,feat_attr,method)
    })
    df["Refseq_idx"]=seq_name
    df.to_csv("Temp_incremental_mutation/sequences_mutated.csv",index=False,header=False,mode="a")


def write_all_incrementally_mutated_seqs(seq_file,seq_colname,feat_attr_file,method):
    """
    Write all mutated sequences of all sequences in seq_file into Temp_incremental_mutation/sequences_mutated.csv
    """
    colnames=["mutated_sequence","Refseq_idx"]

    # write header
    with open("Temp_incremental_mutation/sequences_mutated.csv","w+") as f:
        f.write(",".join(colnames)+"\n")

    seqs=pd.read_csv(seq_file).loc[:,seq_colname]
    attrs=pd.read_csv(feat_attr_file,index_col=0,header=None).values
    for i in range(len(seqs)):
        write_incremental_mutated_sequences_for_one_seq(seqs[i],attrs[i],"Ref"+str(i),method)



def calculate_incremental_mutation(device,model,seq_file,seq_colname,feat_attr_file,key,out_path,method):
    """ Sequence needs to be 600 bp"""
    # make directory for temp result
    if not os.path.isdir("Temp_incremental_mutation/"):
        os.makedirs("Temp_incremental_mutation/")

    # write all mutated sequences
    write_all_incrementally_mutated_seqs(seq_file,seq_colname,feat_attr_file,method)

    # write signal for mutated sequences
    write_predictions(device,model,"Temp_incremental_mutation/sequences_mutated.csv",
                            "mutated_sequence","Temp_incremental_mutation/prediction_mutated.csv",key)

    # find transition point for each sequence signal
    pred=pd.read_csv("Temp_incremental_mutation/prediction_mutated.csv").loc[:,"pred_"+key].values.reshape(-1,600)
    transition_idx=[find_first_two_consecutive_negs(pred[i]) for i in range(pred.shape[0])]
    pd.DataFrame(transition_idx).to_csv(out_path,index=False,header=False)
    shutil.rmtree("Temp_incremental_mutation/")



def incremental_mutations_guided_by_feat_attr(device,model,seq,feat_attr,key,method):
    """
    Mainly used for plotting change of z score
    Performing incremental mutation for single sequence.
    Order mutation locations by decreasing importance.
    Then incrementally mutate each location to "N".
    Finally return signal,z,and pred for all the mutated sequences
    """
    mutated_seqs=incremental_mutated_sequences_for_one_seq(seq,feat_attr,method)
    signal,z,pred=predict_from_seq(model,mutated_seqs,device,key)
    idx=find_first_two_consecutive_negs(pred.squeeze()) 
    return signal,z,pred







#-------------------------
# consecutively mutate to worst SNP
#------------------------
def mutate_worst_SNP_for_one_seq(seq,attr):
    conversion_dict={0:"A",1:"C",2:"G",3:"T"}
    if attr.shape[-1]==2400:
        attr=reshape_attr(attr)
    attrxinp=calculate_attrxinp(seq,attr)
    location=np.argmax(attrxinp)
    attr_at_location=attr.iloc[:,location]
    alt_base=conversion_dict[np.argmin(attr_at_location)]
    mutated_seq=mutate(seq,alt_base,location)
    return mutated_seq



def one_mutation_iteration(device,model,seq_file,num_seqs,seq_colname,num_iteration,target):
    """
    One mutation iteration
    Attribution is calculated from seq_file
    Mutated sequences are appended to existing file
    """
    if num_iteration==0:
        seqs=pd.read_csv(seq_file).loc[:,seq_colname].tolist()
    else:
        seqs_df=pd.read_csv("Temp_worst_SNP/sequences_mutated.csv").iloc[-num_seqs:,:]
        seqs=seqs_df.loc[:,"mutated_sequence"].tolist()
    attr_file_name=f'Temp_worst_SNP/attributions_iteration_{num_iteration}.csv'
    # write attributions in mini-batch manner
    write_attr(device,model,seq_file,seq_colname,"Gradient",attr_file_name,target)
    attrs=pd.read_csv(attr_file_name,index_col=0,header=None)
    # mutate sequences
    mutated_seqs=[mutate_worst_SNP_for_one_seq(seqs[i],attrs.iloc[i,:]) for i in range(len(seqs))]
    df=pd.DataFrame(mutated_seqs,columns=["mutated_sequence"])
    df["mutation_interation"]=num_iteration
    df.to_csv("Temp_worst_SNP/sequences_mutated.csv",header=False,mode="a")


def incremental_worst_mutation(device,model,seq_file,seq_colname,key,iterations,out_path):
    """ Sequence needs to be 600 bp"""
    # make directory for temp result
    if not os.path.isdir("Temp_worst_SNP/"):
        os.makedirs("Temp_worst_SNP/")

    write_predictions(device,model,seq_file,
                           seq_colname,"Temp_worst_SNP/prediction_reference.csv",key,
                           variable_length=True,batch_size=8192)

    num_seqs=pd.read_csv(seq_file).shape[0]

    # write all mutated sequences
    index_dict={"cage_hepg2":0,"cage_k562":1}
    with open("Temp_worst_SNP/sequences_mutated.csv","w+") as f:
            f.write(",".join(["mutated_sequence","mutation_interation"])+"\n")
    for i in range(iterations):
        one_mutation_iteration(device,model,seq_file,num_seqs,seq_colname,i,index_dict[key])

    # write signal for mutated sequences
    write_predictions(device,model,"Temp_worst_SNP/sequences_mutated.csv",
                           "mutated_sequence","Temp_worst_SNP/prediction_mutated.csv",key)


    z_alt=pd.read_csv("Temp_worst_SNP/prediction_mutated.csv").loc[:,f"z_{key}"].values
    z_ref=pd.read_csv("Temp_worst_SNP/prediction_reference.csv").loc[:,f"z_{key}"].values
    zs=np.concatenate((z_ref,z_alt))
    res=zs.reshape(-1,iterations+1,order='F')
    pd.DataFrame(res).to_csv(out_path,index=False)
    shutil.rmtree("Temp_worst_SNP/")

