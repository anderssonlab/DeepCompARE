#!/usr/bin/env python3
import pandas as pd
import os
import pyranges as pr
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor, ablate_motifs
from region_ops import resize_df
from prediction import compute_predictions
from seq_annotators import JasparAnnotator, gnomadAnnotator, phylopAnnotator
#-----------------------
# Functions
#-----------------------

# TODO: refactor, and add region info, use get_motif_isa(), make gnomadAnnotator smaller, and remove it and gc() when done

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
                        phylop241way_annotator,
                        phylop447way_annotator,
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
    motif_df = get_phylop(motif_df,"241way",phylop241way_annotator)
    motif_df = get_phylop(motif_df,"447way",phylop447way_annotator)
    motif_df = motif_df.round(3)

    # add sequence information
    motif_df["seq_idx"]= f"Seq{idx}"
    motif_df["region"]= f"{region[0]}:{region[1]}-{region[2]}"
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
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{args.file_name}.bed").df
    # rename columns
    df_regions.columns=["chromosome","start","end","name"]
    df_regions=resize_df(df_regions,600)
    df_regions["end"]=df_regions["end"]-1
    
    # for debug purposes
    df_regions=df_regions.iloc[:10,:]
    
    # load sequence extractor
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")

    # load jaspar annotator
    if "hepg2" in args.file_name:
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
    if "k562" in args.file_name:
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                     score_thresh=500,
                                     rna_file=rna_file)
    
    df_motif=jaspar_annotator.annotate(df_regions)
    
    gnomad_annotator=gnomadAnnotator()
    phylop241way_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP241way.bw")
    phylop447way_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP447way.bw")
    
    
    
    
    logger.info(f"Finished annotating {args.file_name}.")

