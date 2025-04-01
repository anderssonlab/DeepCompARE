#!/usr/bin/env python3
import pandas as pd
import torch
import os
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator, gnomadAnnotator, phylopAnnotator, baseImpAnnotator
from in_silico_mutagenesis import get_motif_isa

#-----------------------
# Functions
#-----------------------

THRESHOLD=500

def main(file_name,device):
    # Load data and tools
    df_regions=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/{file_name}.tsv",sep=" ").iloc[:,0:3]
    df_regions.columns=["chromosome","start","end"]
    df_regions["end"]=df_regions["end"]-2
    # select chromosomes
    df_regions=df_regions[df_regions["chromosome"].apply(lambda x: x in [f"chr{i}" for i in list(range(1,23))+["X","Y"]])].reset_index(drop=True)
    
    # Step 1: Annotate motif info
    if "hepg2" in file_name.lower():
        chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed"
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
    if "k562" in file_name:
        chip_file= "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed"
        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                        score_thresh=THRESHOLD,
                                        chip_file=chip_file,
                                        rna_file=rna_file)
    jaspar_annotator.annotate(df_regions,outpath=f"{file_name}_temp1.csv")
    logger.info("Done with Step 1: motif extraction")

    # step 2: read back motif info, get isa
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    df_motif=pd.read_csv(f"{file_name}_temp1.csv")
    df_motif=get_motif_isa(seq_extractor,df_motif,device=torch.device(f"cuda:{device}"))
    logger.info(f"Done with gpu for file {file_name}")
    df_motif.to_csv(f"{file_name}_temp2.csv",index=False)
    
    # step 3: annotate with conservation and feature importance
    df_motif=pd.read_csv(f"{file_name}_temp2.csv")
    phylop_annotator=phylopAnnotator(f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP241way.bw")
    gradxinp_annotator=baseImpAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/gradxinp_{file_name}.csv")
    ism_annotator=baseImpAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/ism_{file_name}.csv")
    isa_annotator=baseImpAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/isa_{file_name}.csv")
    
    for chr_num in list(range(1,23))+["X","Y"]: 
        logger.info(f"Working on chromosome {chr_num}")
        df_chr=df_motif[df_motif["chromosome"]==f"chr{chr_num}"].reset_index(drop=True)
        if df_chr.shape[0]==0:
            continue
        gnomad_annotator=gnomadAnnotator(chr_num)
        df_chr = gnomad_annotator.annotate(df_chr)
        df_chr["phylop_241way"] = df_chr.apply(lambda row: phylop_annotator.annotate((row['chromosome'],row['start'],row['end'])), axis=1)
        for track_num in range(8):
            logger.info(f"Working on track {track_num}")
            df_chr[f"gradxinp_{track_num}"]=df_chr.apply(lambda row: gradxinp_annotator.annotate(row['region'],track_num,row['start_rel'],row['end_rel']), axis=1)
            df_chr[f"ism_{track_num}"]=df_chr.apply(lambda row: ism_annotator.annotate(row['region'],track_num,row['start_rel'],row['end_rel']), axis=1)
            df_chr[f"isa_{track_num}"]=df_chr.apply(lambda row: isa_annotator.annotate(row['region'],track_num,row['start_rel'],row['end_rel']), axis=1)
        # if the file doesn't exist, write header
        out_name=f"motif_info_thresh_{THRESHOLD}_{file_name}.csv"
        if not os.path.isfile(out_name):
            df_chr.to_csv(out_name,mode="w",header=True)
        else:
            df_chr.to_csv(out_name,mode="a",header=False)
    
    # remove temp file
    # os.remove(f"{file_name}_temp1.csv")
    # os.remove(f"{file_name}_temp2.csv")



#-----------------------
# Analysis
#-----------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    parser.add_argument('--device', type=str)
    
    args=parser.parse_args()
    main(args.file_name,args.device)
    logger.info(f"Finished annotating {args.file_name}.")

