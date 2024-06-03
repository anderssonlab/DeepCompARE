import pandas as pd
import pyranges as pr
import sys
from loguru import logger
import argparse
import os

sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor
from variant_interpreter import predict_vcf
from seq_annotators import JasparAnnotator, ReMapAnnotator, read_maf_file



def subset_maf_by_range(range_bed_file, gr_maf):
    gr_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{range_bed_file}.bed")
    df_maf = gr_maf.join(gr_regions).as_df()
    df_maf.drop(["Start_b","End_b"],axis=1,inplace=True)
    df_maf.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']
    return df_maf



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Correlation analysis")
    parser.add_argument("--range_bed_file",type=str,help="bed file for range")
    parser.add_argument("--gpu",type=str,help="gpu number")
    args=parser.parse_args()
    
    # get remap annotator
    jaspar_annotator=JasparAnnotator()
    if "k562" in args.range_bed_file:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")
    if "hepg2" in args.range_bed_file:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_Hep-G2.bed")
    
    # process one chromosome {args.chrom} on one range (promoter/enhancer) {args.range_bed_file} 
    for i in list(range(1,23))+["X","Y"]:
        df_maf=read_maf_file(i) 
        gr_maf = pr.PyRanges(df_maf)
        df_maf=subset_maf_by_range(args.range_bed_file,gr_maf)
        
        # get variant effect prediction using deepCompare
        seq_extractor = SeqExtractor()
        df_pred=pd.DataFrame(predict_vcf(df_maf,seq_extractor,args.gpu))
        df_res=pd.concat([df_maf,df_pred],axis=1)
        df_res.to_csv(f"maf_with_effect_size_{args.range_bed_file}.csv",mode="a",header=False,index=False)
        
    logger.info(f"Done with {args.range_bed_file}")
