"""
This module takes Pd5_motif_info/motif_info_XXX as input, 
calculate the mutational constraints of each TF from column "af"
"""


import pandas as pd
import numpy as np



def get_mutation_dict():
    df=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_constraints/Supplementary_Data_1.tsv", sep="\t")
    # group by cotext, ref, alt, sum up possible_variants, observed_variants
    df=df.groupby(["context","ref","alt"]).agg({"possible_variants":"sum","observed_variants":"sum"}).reset_index()
    # group by context, avg the possible_variants, sum the observed_variants    
    df=df.groupby("context").agg({"possible_variants":"mean","observed_variants":"sum"}).reset_index()
    df["mutation_rate"]=df["observed_variants"]/df["possible_variants"]
    # add reverse complement
    df_rc=df.copy()
    df_rc["context"]=df_rc["context"].apply(lambda x: x[::-1].translate(str.maketrans("ACGT","TGCA")))
    # merge the two dfs
    df=pd.concat([df,df_rc], axis=0)
    # drop columns possible_variants, observed_variants
    df.drop(["possible_variants","observed_variants"], axis=1, inplace=True)
    # convert to dict
    mutation_dict=df.set_index("context").to_dict()["mutation_rate"]
    return mutation_dict



def num_expected_mutations(motif, mutation_dict):
    # step 1: get all trinucleotides in the motif
    trinucleotides = [motif[i:i+3] for i in range(len(motif)-2)]
    expected_mutations = sum([mutation_dict[tri] for tri in trinucleotides])
    return expected_mutations


def calc_gnocchi(df):
    df['oe'] = df['num_rare_variants']/df['num_expected_variants']
    df['chisq'] = (df['num_rare_variants']-df['num_expected_variants'])**2 / df['num_expected_variants']
    df['z'] = np.where(df['oe'] >= 1., (-1)*np.sqrt(df['chisq']), np.sqrt(df['chisq']))
    # convert larger z to 10
    df['z'] = np.where(df['z'] > 10, 10, df['z'])
    df['z'] = np.where(df['z'] < -10, -10, df['z'])
    df.reset_index(drop=True, inplace=True)
    return df




def calc_constraint(df_orig,seq_extractor):
    df=df_orig.copy()
    mutation_dict=get_mutation_dict()
    df["motif"]=df.apply(lambda row: seq_extractor.get_seq((row["chromosome"],row["start"],row["end"])), axis=1)
    df["num_rare_variants"] = df["gnomad_af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["num_expected_variants"]=df["motif"].apply(lambda x: num_expected_mutations(x, mutation_dict))
    df=df.groupby("protein").agg({"num_expected_variants":"sum",
                                  "num_rare_variants":"sum"}).reset_index()
    df=calc_gnocchi(df)
    return df

