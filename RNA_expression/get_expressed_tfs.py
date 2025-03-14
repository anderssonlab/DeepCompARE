import pandas as pd




thresh=0.5


def map_name(df):
    df_name_mapper=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv")
    df_name_mapper=df_name_mapper[['Ensembl ID', 'HGNC symbol']].copy()
    df=pd.merge(df,df_name_mapper,left_on='gene_id',right_on='Ensembl ID',how='inner')
    return df


def get_expressed_proteins(file_path):
    df_expr=pd.read_csv(file_path,sep='\t')
    df_expr["gene_id"]=df_expr["gene_id"].str.split(".").str[0]
    df_expr=map_name(df_expr)
    # select the HGNC symbol and TPM>thresh
    df_expr=df_expr[df_expr['TPM']>thresh].reset_index(drop=True)
    return df_expr['HGNC symbol'].tolist()



def write_expressed_tfs(rep_path1,
                        rep_path2,
                        output_path,
                        ):
    proteins1_rna=get_expressed_proteins(rep_path1)
    proteins2_rna=get_expressed_proteins(rep_path2)
    # get the proteins expressed in both replicates
    proteins_rna=list(set(proteins1_rna).intersection(set(proteins2_rna)))
    protein_list=proteins_rna 
    pd.DataFrame(protein_list).to_csv(output_path,header=None,index=None)




write_expressed_tfs("Raw_data/ENCFF928NYA_rep1_RNASeq_K562.tsv",
                    "Raw_data/ENCFF003XKT_rep2_RNASeq_K562.tsv",
                    "expressed_tf_list_k562.tsv")
write_expressed_tfs("Raw_data/ENCFF103FSL_rep1_RNASeq_HepG2.tsv",
                    "Raw_data/ENCFF692QVJ_rep2_RNASeq_HepG2.tsv",
                    "expressed_tf_list_hepg2.tsv")



# Final decision: 
# threshold=0.5
# no chip-seq
# don't subset by htf because it doesn't have SMAD2