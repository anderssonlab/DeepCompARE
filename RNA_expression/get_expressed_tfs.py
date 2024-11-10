import pandas as pd


def map_name(df):
    df_name_mapper=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv")
    df_name_mapper=df_name_mapper[['Ensembl ID', 'HGNC symbol']].copy()
    df=pd.merge(df,df_name_mapper,left_on='gene_id',right_on='Ensembl ID',how='inner')
    return df


def write_expressed_tfs(file_path, chip_path, output_path):
    df_expr=pd.read_csv(file_path,sep='\t')
    df_expr["gene_id"]=df_expr["gene_id"].str.split(".").str[0]
    df_expr=map_name(df_expr)
    # select the HGNC symbol and TPM>1
    df_expr=df_expr[df_expr['TPM']>1].reset_index(drop=True)
    # output the HGNC symbol
    protein_names=df_expr['HGNC symbol'].tolist()
    chip_tfs=pd.read_csv(chip_path,header=None).iloc[:,0].tolist()
    protein_list=list(set(protein_names+chip_tfs))
    pd.DataFrame(protein_list).to_csv(output_path,header=None,index=None)



write_expressed_tfs("Raw_data/ENCFF928NYA_RNASeq_K562.tsv",
                    "/isdata/alab/people/pcr980/Resource/ReMap2022/TFs_K562.txt",
                    "expressed_tf_list_k562.tsv")
write_expressed_tfs("Raw_data/ENCFF103FSL_RNASeq_HepG2.tsv",
                    "/isdata/alab/people/pcr980/Resource/ReMap2022/TFs_HepG2.txt",
                    "expressed_tf_list_hepg2.tsv")

