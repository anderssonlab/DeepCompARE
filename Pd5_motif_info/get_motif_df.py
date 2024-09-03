#!/usr/bin/env python3
import pandas as pd
import os
import pyranges as pr
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor, ablate_motifs
from prediction import compute_predictions
from seq_annotators import JasparAnnotator, BindingEvidenceAnnotator, gnomadAnnotator, phylopAnnotator
#-----------------------
# Functions
#-----------------------


def get_predictions(seq_col,device,prefix):
    pred=compute_predictions(seq_col,device=device)
    df_pred=pd.DataFrame(pred,columns=[f"pred_{prefix}_track{i}" for i in range(16)])
    return df_pred




def get_isms(seq_extractor,motif_df,region,device):
    sequence=seq_extractor.get_seq(region)
    motif_df["seq_orig"]= sequence
    motif_df['seq_mut'] = motif_df.apply(lambda row: ablate_motifs(sequence, row['start_rel'], row['end_rel']), axis=1)
    df_pred_orig=get_predictions(motif_df["seq_orig"],device=device,prefix="orig").round(3)
    df_pred_mut=get_predictions(motif_df["seq_mut"],device=device,prefix="mut").round(3)
    motif_df=pd.concat([motif_df,df_pred_orig,df_pred_mut],axis=1)
    for i in range(16):
        motif_df[f"ism_track{i}"]=motif_df[f"pred_orig_track{i}"]-motif_df[f"pred_mut_track{i}"]
        motif_df[f"ism_track{i}"]=motif_df[f"ism_track{i}"].round(3)
        motif_df.drop(columns=[f"pred_mut_track{i}"],inplace=True)
    motif_df.drop(columns=["seq_orig","seq_mut"],inplace=True)
    return motif_df



def get_phylop(motif_df,prefix, phylop_annotator):
    motif_df[f"{prefix}"] = motif_df.apply(lambda row: phylop_annotator.annotate((row['chromosome'],row['start'],row['end'])), axis=1)
    return motif_df



# TODO: change threshold here
def annotate_one_region(idx,
                        region,
                        out_path,
                        seq_extractor,
                        jaspar_annotator, 
                        be_annotator,
                        gnomad_annotator,
                        phylop100way_annotator,
                        phylop241way_annotator,
                        phylop447way_annotator,
                        phylop447wayLRT_annotator, 
                        device,
                        score_thresh=500):
    motif_df=jaspar_annotator.annotate(region)
    motif_df=motif_df[motif_df["score"]>score_thresh].reset_index(drop=True)
    if motif_df.shape[0]==0:
        return
    # add relative start and end
    logger.info(f"Annotating {motif_df.shape[0]} motifs in region {idx}.")
    motif_df["start_rel"]=motif_df["start"]-region[1]
    motif_df["end_rel"]=motif_df["end"]-region[1]
    # add whether it has binding evidence
    motif_df = be_annotator.annotate(motif_df,region)
    # add ism
    motif_df = get_isms(seq_extractor,motif_df,region,device)
    # add gnomad info
    motif_df["af"] = motif_df.apply(lambda row: gnomad_annotator.annotate((row['chromosome'],row['start'],row['end'])), axis=1)
    # add phylop info
    motif_df = get_phylop(motif_df,"100way",phylop100way_annotator)
    motif_df = get_phylop(motif_df,"241way",phylop241way_annotator)
    motif_df = get_phylop(motif_df,"447way",phylop447way_annotator)
    motif_df = get_phylop(motif_df,"447wayLRT",phylop447wayLRT_annotator)
    motif_df = motif_df.round(3)

    # add sequence information
    motif_df["seq_idx"]= f"Seq{idx}"
    # write to csv
    if not os.path.exists(out_path):
        motif_df.to_csv(out_path,index=False,mode="w",header=True)
    else:
         motif_df.to_csv(out_path,index=False,mode="a",header=False)



#-----------------------
# Analysis
#-----------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    parser.add_argument('--device', type=str)
    
    args=parser.parse_args()
    
    # Load data and tools
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/promoters_hepg2.bed").df
    #df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/{args.file_name}.bed").df
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    
    # load tools
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",by="contained")
    gnomad_annotator=gnomadAnnotator()
    phylop100way_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP100way.bw")
    phylop241way_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP241way.bw")
    phylop447way_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP447way.bw")
    phylop447wayLRT_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP447wayLRT.bw")
    
    if "k562" in args.file_name.lower():
        be_annotator = BindingEvidenceAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_K-562.bed",mode="chip")
    if "hepg2" in args.file_name.lower():
        be_annotator = BindingEvidenceAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_Hep-G2.bed",mode="chip")

    for idx in range(df_regions.shape[0]):
        if idx%1000==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,0:3].tolist()
        # TODO: change file name relative to "thresh" here
        annotate_one_region(idx,
                            region,
                            f"motif_info_thresh_500_{args.file_name}.csv",
                            seq_extractor,
                            jaspar_annotator,
                            be_annotator,
                            gnomad_annotator,
                            phylop100way_annotator,
                            phylop241way_annotator,
                            phylop447way_annotator,
                            phylop447wayLRT_annotator,
                            args.device)
    logger.info(f"Finished annotating {args.file_name}.")

