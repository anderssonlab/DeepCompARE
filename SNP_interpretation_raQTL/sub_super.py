import pandas as pd
import numpy as np
import pyranges as pr



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_annotators import JasparAnnotator, ReMapAnnotator

cell_type="k562"
df_processed=pd.read_csv(f"Raw_data/{cell_type}_sure_raQTL.csv")
df_processed.drop("ID",axis=1,inplace=True)
# make sure SNP_ID is unique
df_processed=df_processed.drop_duplicates(subset="SNP_ID").reset_index(drop=True)
assert df_processed.SNP_ID.nunique()==df_processed.shape[0]
df_processed.set_index("SNP_ID",inplace=True)


df_raw=pd.read_csv(f"Raw_data/{cell_type}.sign.id.LP190708.txt",sep="\t")
assert np.all(df_processed.index.isin(df_raw.SNP_ID))
assert df_raw.SNP_ID.nunique()==df_raw.shape[0]
# set SNP_ID as index
df_raw.set_index("SNP_ID",inplace=True)

# merge
df_processed=df_processed.join(df_raw,how="inner",rsuffix="_raw")
np.all(df_processed.ref==df_processed.ref_raw)
df_processed=df_processed[["chr","SNPabspos","ref","alt","log2FC"]]

# convert df_processed to pyranges
df_processed["Start"]=df_processed.SNPabspos-1
df_processed["End"]=df_processed.SNPabspos
gr_qtl=pr.PyRanges(df_processed.rename(columns={"chr":"Chromosome"}))


promoters=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/promoters_{cell_type}.bed")
df_qtl_promoter=gr_qtl.overlap(promoters).df
# re-order columns
df_qtl_promoter=df_qtl_promoter[["Chromosome","Start","End","ref","alt","log2FC"]]

enhancers=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/enhancers_{cell_type}.bed")
df_qtl_enhancer=gr_qtl.overlap(enhancers).df
# re-order columns
df_qtl_enhancer=df_qtl_enhancer[["Chromosome","Start","End","ref","alt","log2FC"]]


# annotate withh Jaspar and ReMap
jaspar_annotator=JasparAnnotator(jaspar_path="/binf-isilon/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg19.bb")
remap_annotator=ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")

for i in range(df_qtl_enhancer.shape[0]):
    region=df_qtl_enhancer.iloc[i,0:3]
    motifs=jaspar_annotator.annotate(region)
    print(motifs) # all empty, end of search