"""
       Promoter Common Enhancer
HepG2
Common
K562

Column promoter: TFs appear in promoters but not enhancers, regardless of cell type

"""


import pandas as pd

df=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_targets_and_properties/tf_individual_effect_by_file.csv')
# use only "dstat_ism_*"
columns=df.columns[df.columns.str.contains('dstat_ism_')].to_list()
df=df[columns+["protein","dataset"]]
df.corr() # very high, so avg or PC1 is reasonable
df["assay_avg"]=df[columns].mean(axis=1)
df=df[["protein","dataset","assay_avg"]]


promoter_tfs=df[df.dataset.str.contains('promoter') & df.assay_avg>0]["protein"].unique()
enhancer_tfs=df[df.dataset.str.contains('enhancer') & df.assay_avg>0]["protein"].unique()
common_pe_tfs=set(promoter_tfs).intersection(set(enhancer_tfs))



