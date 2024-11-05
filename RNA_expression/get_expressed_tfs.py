import pandas as pd


def map_name(df):
    df_name_mapper=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_transcription_factors]/DatabaseExtract_v_1.01.csv")
    df_name_mapper=df_name_mapper[['Ensembl ID', 'HGNC symbol']].copy()
    df=pd.merge(df,df_name_mapper,left_on='gene_id',right_on='Ensembl ID',how='inner')
    return df


def write_expressed_tfs(file_path, output_path):
    df_expr=pd.read_csv(file_path,sep='\t')
    df_expr["gene_id"]=df_expr["gene_id"].str.split(".").str[0]
    df_expr=map_name(df_expr)
    # select the HGNC symbol and TPM>1
    df_expr=df_expr[df_expr['TPM']>1].reset_index(drop=True)
    # output the HGNC symbol
    df_expr['HGNC symbol'].to_csv(output_path,index=False,header=False)



write_expressed_tfs("Raw_data/ENCFF928NYA_RNASeq_K562.tsv","expressed_tfs_k562.txt")
write_expressed_tfs("Raw_data/ENCFF103FSL_RNASeq_HepG2.tsv","expressed_tfs_hepg2.txt")

