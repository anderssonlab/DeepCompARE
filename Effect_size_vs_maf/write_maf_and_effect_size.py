import pandas as pd
import pyranges as pr
import sys
from loguru import logger
import argparse

sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor
from variant_interpreter import predict_vcf
from motif_annotation import JasparAnnotator, ReMapAnnotator



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
    gr_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{range_bed_file}.bed")
    df_maf = gr_maf.join(gr_regions).as_df()
    df_maf.drop(["Start_b","End_b"],axis=1,inplace=True)
    df_maf.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']
    return df_maf

    
def get_top_score_tf(df_motif):
    """
    Args:
    df_motif: dataframe with motif annotation
    
    Return:
    tf: TF with highest score
        prioritize TFBSs with ChIP evidence with highest score
        if no TFBS has ChIP evidence, then return the one with highest score
    """

    df_chip_true=df_motif[df_motif["chip_evidence"]==True]
    if df_chip_true.shape[0]>0:
        df_chip_true=df_chip_true.sort_values("score",ascending=False).reset_index()
        tf=df_chip_true.protein[0]
        score=df_chip_true.score[0]
    else:
        df_motif=df_motif.sort_values("score",ascending=False).reset_index()
        tf=df_motif.protein[0]
        score=df_motif.score[0]
    return tf,score


def annotate_one_row(row,jaspar_annotator,remap_annotator):
    region=(row[0],row[1],row[2])
    df_motif=jaspar_annotator.annotate(region,by="1bp")
    if df_motif.shape[0]==0:
        return None, None
    df_motif = remap_annotator.annotate(df_motif,region)
    return get_top_score_tf(df_motif)



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Correlation analysis")
    parser.add_argument("--range_bed_file",type=str,help="bed file for range")
    parser.add_argument("--gpu",type=str,help="gpu number")
    parser.add_argument("--chrom",type=str,help="chromosome number")
    args=parser.parse_args()
    
    # get remap annotator
    jaspar_annotator=JasparAnnotator()
    if "k562" in args.range_bed_file:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")
    if "hepg2" in args.range_bed_file:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_Hep-G2.bed")
    
    # process one chromosome {args.chrom} on one range (promoter/enhancer) {args.range_bed_file} 
    df_maf=read_maf_file(args.chrom) 
    gr_maf = pr.PyRanges(df_maf)
    df_maf=subset_maf_by_range(args.range_bed_file,gr_maf)
    
    tfbs_list=[]
    motif_score_list=[]
    for i,row in df_maf.iterrows():
        if i%1000==0:
            logger.info(f"Processing {i}th row")
        tfbs,motif_score=annotate_one_row(row,jaspar_annotator,remap_annotator)
        tfbs_list.append(tfbs)
        motif_score_list.append(motif_score)
    
    df_maf["TFBS"]=tfbs_list
    df_maf["motif_score"]=motif_score_list
    
    # get variant effect prediction using deepCompare
    seq_extractor = SeqExtractor()
    df_pred=pd.DataFrame(predict_vcf(df_maf,seq_extractor,args.gpu))
    df_res=pd.concat([df_maf[["AF","TFBS","motif_score"]],df_pred],axis=1)
    df_res.to_csv(f"chrom{args.chrom}.csv",mode="a",header=False,index=False)
        
    logger.info(f"Done with {args.range_bed_file}")
