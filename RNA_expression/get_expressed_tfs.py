import pandas as pd




thresh=1

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



def write_expressed_tfs(rep_path1,rep_path2, chip_path, output_path):
    proteins1=get_expressed_proteins(rep_path1)
    proteins2=get_expressed_proteins(rep_path2)
    # get the proteins that are expressed in both replicates
    proteins=list(set(proteins1).intersection(set(proteins2)))
    # output the HGNC symbol
    # chip_tfs=pd.read_csv(chip_path,header=None).iloc[:,0].tolist()
    # protein_list=list(set(proteins+chip_tfs))
    protein_list=proteins
    pd.DataFrame(protein_list).to_csv(output_path,header=None,index=None)
    return protein_list




write_expressed_tfs("Raw_data/ENCFF928NYA_rep1_RNASeq_K562.tsv",
                    "Raw_data/ENCFF003XKT_rep2_RNASeq_K562.tsv",
                    "/isdata/alab/people/pcr980/Resource/ReMap2022/TFs_K562.txt",
                    f"expressed_tf_list_k562_{thresh}.tsv")
write_expressed_tfs("Raw_data/ENCFF103FSL_rep1_RNASeq_HepG2.tsv",
                    "Raw_data/ENCFF692QVJ_rep2_RNASeq_HepG2.tsv",
                    "/isdata/alab/people/pcr980/Resource/ReMap2022/TFs_HepG2.txt",
                    f"expressed_tf_list_hepg2_{thresh}.tsv")

