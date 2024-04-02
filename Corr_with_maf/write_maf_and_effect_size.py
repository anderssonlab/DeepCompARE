import pandas as pd
import numpy as np
import pyranges as pr
import sys
from loguru import logger
import argparse

sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor
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


def subset_maf_by_range(range_bed_file, gr_maf):
    # get range
    gr_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{range_bed_file}.bed")
    df_maf = gr_maf.join(gr_regions).as_df()
    df_maf.drop(["Start_b","End_b"],axis=1,inplace=True)
    df_maf.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']
    return df_maf
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Correlation analysis")
    parser.add_argument("--range_bed_file",type=str,help="bed file for range")
    args = parser.parse_args()
    
    for i in list(range(1,23))+["X","Y"]:
        logger.info(f"Processing chromosome {i}")
        df_maf=read_maf_file(i)
        gr_maf = pr.PyRanges(df_maf)
        df_maf=subset_maf_by_range(args.range_bed_file,gr_maf)

        seq_extractor = SeqExtractor()
        df_pred=pd.DataFrame(predict_vcf(df_maf,seq_extractor))
        df_res=pd.concat([df_maf[["AF"]],df_pred],axis=1)
        df_res.to_csv(f"maf_with_effect_size_{args.range_bed_file}.csv",mode="a",header=False,index=False)
        
    logger.info(f"Done with {args.range_bed_file}")