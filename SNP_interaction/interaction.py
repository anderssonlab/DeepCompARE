from Bio import SeqIO
import pandas as pd
import sys
import os
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
import numpy as np
from write_prediction import compute_predictions,write_predictions

#---------------------
# Task1: get sequences from hg38
#---------------------

# Define genomic regions
# regions = [
#     ("chr10", 132396935, 132397534),
#     ("chr10", 8054655, 8055254),
#     ("chr15", 37098712, 37099311),
#     ("chr12", 132828699, 132829298),
#     ("chr13", 114234732, 114235331),
#     ("chr17", 36544867, 36545466)
# ]


# def get_sequence(chrom, start, end):
#     sequence = hg38_genome[chrom].seq[start-1:end]
#     return str(sequence).upper()

# with open("/isdata/alab/people/pcr980/Hg38/hg38.fa", "r") as file:
#     hg38_genome = SeqIO.to_dict(SeqIO.parse(file,"fasta"))
#     seqs=[get_sequence(*region) for region in regions]
    
# seqs_df=pd.DataFrame(seqs)
# seqs_df.columns=["sequence"]
# seqs_df.to_csv("/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/seqs_ref.csv",index=False)




#---------------------
# Task2: 
# create table with:
# 1st mutation location, 2nd mutation location, 
# 1st location reference base, 2nd location reference base,
# 1st mutation base, 2nd mutation base, 
# sequence with 1st mutation, sequence with 2nd mutation, sequence with both mutations
# calculate prediction for each mutated sequence, compare to reference
#---------------------

def get_alt_bases(base_ref):
    alt_bases=["A","C","G","T"]
    alt_bases.remove(base_ref)
    return alt_bases

def mutate_single_base(seq,alt_base,location):
    return seq[:location]+alt_base+seq[location+1:]

def mutate_2_bases(seq,alt_base1,alt_base2,location1,location2):
    return seq[:location1]+alt_base1+seq[location1+1:location2]+alt_base2+seq[location2+1:]

def test_mutate_2_bases():
    seq="AAAAAAAA"
    alt_base1="G"
    alt_base2="T"
    location1=1
    location2=3
    print(mutate_2_bases(seq,alt_base1,alt_base2,location1,location2)=="AGATAAAA")
    return 



def process_one_seq(seq,seq_name):
    # create empty dataframe
    df_info=pd.DataFrame(columns=["mut1_loc","mut2_loc","mut1_ref","mut2_ref","mut1_alt","mut2_alt","seq_mut1","seq_mut2","seq_mut1_mut2"])   
    
    # Generate all possible mutations      
    mut_locs = range(600)
    all_combinations = [(i, j) for i in mut_locs for j in mut_locs if j > i]

    # Generate mutations
    mutations = [(mut1_loc, mut1_alt, mut2_loc, mut2_alt)
                for mut1_loc, mut2_loc in all_combinations
                for mut1_alt in get_alt_bases(seq[mut1_loc])
                for mut2_alt in get_alt_bases(seq[mut2_loc])]

    # Convert to DataFrame
    mut_cols = ["mut1_loc", "mut1_alt", "mut2_loc", "mut2_alt"]
    df_info = pd.DataFrame(mutations, columns=mut_cols)

    # Write mutated sequences in vectorized form
    df_info["seq_mut1"] = np.vectorize(mutate_single_base)(seq, df_info["mut1_alt"], df_info["mut1_loc"])
    df_info["seq_mut2"] = np.vectorize(mutate_single_base)(seq, df_info["mut2_alt"], df_info["mut2_loc"])
    df_info["seq_mut1_mut2"] = np.vectorize(mutate_2_bases)(seq, df_info["mut1_alt"], df_info["mut2_alt"], df_info["mut1_loc"], df_info["mut2_loc"])

    # Add reference bases
    df_info["mut1_ref"] = [seq[i] for i in df_info["mut1_loc"]]
    df_info["mut2_ref"] = [seq[i] for i in df_info["mut2_loc"]]

    # Reorder columns to match original structure
    df_info = df_info[["mut1_loc", "mut2_loc", "mut1_ref", "mut2_ref", "mut1_alt", "mut2_alt", "seq_mut1", "seq_mut2", "seq_mut1_mut2"]]       
    
    # write to csv
    df_info.to_csv(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/df_info_{seq_name}.csv",index=False)
   
    # Calculate prediction for each mutated sequence
    write_predictions(data_path=f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/df_info_{seq_name}.csv",
                    seq_colname="seq_mut1",
                    out_path=f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut1_{seq_name}.csv",)

    write_predictions(data_path=f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/df_info_{seq_name}.csv",
                    seq_colname="seq_mut2",
                    out_path=f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut2_{seq_name}.csv",)

    write_predictions(data_path=f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/df_info_{seq_name}.csv",
                    seq_colname="seq_mut1_mut2",
                    out_path=f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut1_mut2_{seq_name}.csv",)
    
    # For now only use track 1
    pred_mut1=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut1_{seq_name}.csv",header=None).values[:,1]
    pred_mut2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut2_{seq_name}.csv",header=None).values[:,1]
    pred_mut1_mut2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut1_mut2_{seq_name}.csv",header=None).values[:,1]
    
    # calculate reference
    pred_ref=compute_predictions(seq)[0][1]
    
    # add predictions to df_info
    df_info["pred_mut1"]=pred_mut1-pred_ref
    df_info["pred_mut2"]=pred_mut2-pred_ref
    df_info["pred_mut1_mut2"]=pred_mut1_mut2-pred_ref
    
    # write final df_info
    # remove sequence info in df_info to save space
    df_info=df_info.drop(columns=["seq_mut1","seq_mut2","seq_mut1_mut2"])
    df_info.to_csv(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/df_info_{seq_name}.csv",index=False)
    
    #remove intermediate files
    os.remove(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut1_{seq_name}.csv")
    os.remove(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut2_{seq_name}.csv")
    os.remove(f"/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/pred_mut1_mut2_{seq_name}.csv")

seqs_ref=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/SNP_interaction/seqs_ref.csv")

for i,seq in enumerate(seqs_ref["sequence"]):
    process_one_seq(seq,"seq"+str(i))