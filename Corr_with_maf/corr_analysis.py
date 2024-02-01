import pandas as pd
import pyranges as pr
import sys

sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import SeqExtractor
from variant_interpreter import predict_vcf



def read_maf_file(chrom_num):
    df_maf=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Gnomad_vcf/Pd1_MAF/chrom{chrom_num}.csv")
    df_maf=df_maf[df_maf["REF"].str.len()==1]
    df_maf=df_maf[df_maf["ALT"].str.len()==1]
    df_maf=df_maf.reset_index(drop=True)
    df_maf["End"]=df_maf["POS"]+1
    df_maf.columns=["Chromosome","Start","ID","REF","ALT","AF","End"]
    df_maf=df_maf[["Chromosome","Start","End","ID","REF","ALT","AF"]]
    return df_maf

def subset_maf_by_range(gr_maf):
    # get range
    df_regions=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/resize_600bp_CAGE_K562.bed",sep="\t",header=None)
    df_regions=df_regions.iloc[:,0:3]
    df_regions.columns=["Chromosome","Start","End"]
    gr_regions = pr.PyRanges(df_regions)
    df_maf = gr_maf.join(gr_regions).as_df()
    df_maf.drop(["Start_b","End_b"],axis=1,inplace=True)
    df_maf.columns=["chromosome","start","end","ID","REF","ALT","AF"]
    return df_maf
    

if __name__=="__main__":
    for i in list(range(1,23))+["X","Y"]:
        df_maf=read_maf_file(i)
        gr_maf = pr.PyRanges(df_maf)
        df_maf=subset_maf_by_range(gr_maf)

        seq_extractor = SeqExtractor()
        df_pred=pd.DataFrame(predict_vcf(df_maf,seq_extractor))
        df_maf=pd.concat([df_maf[["AF"]],df_pred],axis=1)
        df_maf.to_csv(f"/isdata/alab/people/pcr980/DeepCompare/Corr_with_maf/res.csv",mode="a",header=False,index=False)
        
