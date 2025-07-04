import pandas as pd
import numpy as np
import pyranges as pr
import argparse
import sys
import os


sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from seq_annotators import JasparAnnotator
from loguru import logger
from tf_cooperativity import get_permutations,ablate_motifs,get_predictions
#-----------------------
# Functions
#-----------------------

def extract_motif_and_null_info(jaspar_annotator,region):
    """
    Given a genomic region, extract the motif information.
    Args:
        region: A tuple of (chromosome, start, end).
    Returns:
        A data frame of location, TF, Jaspar score and strand.
    """
    df_motif=jaspar_annotator.annotate(region)
    if df_motif.shape[0]==0:
        return df_motif,None
    
    df_motif=df_motif[df_motif["rna_evidence"]].reset_index(drop=True)
    df_motif.drop(columns=['rna_evidence'],inplace=True)
    # calculate null regions: regions non-overlapping putative TFBSs 
    gr_motif=pr.PyRanges(df_motif.loc[:,['chromosome','start','end']].rename(columns={"chromosome":"Chromosome","start":"Start","end":"End"}))
    gr_region= pr.PyRanges(chromosomes=[region[0]], starts=[region[1]], ends=[region[2]])
    gr_null=gr_region.subtract(gr_motif)
    gr_null.End = gr_null.End - 1
    gr_null.Start = gr_null.Start + 1
    gr_null=gr_null[gr_null.lengths()>9]
    
    df_motif["start_rel"]=df_motif["start"]-region[1]
    df_motif["end_rel"]=df_motif["end"]-region[1]
    df_motif["motif_length"]=df_motif["end_rel"]-df_motif["start_rel"]
    return df_motif,gr_null



def get_isas(df):
    for i in range(8):
        df[f"isa1_null_track{i}"]=df[f"pred_orig_track{i}"]-df[f"pred_mut1_null_track{i}"]
        df[f"isa2_track{i}"]=df[f"pred_orig_track{i}"]-df[f"pred_mut2_track{i}"]
        df[f"isa_both_null_track{i}"]=df[f"pred_orig_track{i}"]-df[f"pred_mut_both_null_track{i}"]
    return df

def write_mutation_res_for_one_region(seq_extractor,
                                      jaspar_annotator,
                                      region,
                                      out_path,
                                      region_idx,
                                      device):
    # get motifs and null regions
    df_motif,gr_null=extract_motif_and_null_info(jaspar_annotator,region)
    if df_motif.shape[0]==0:
        return
    
    # get indices of all possible motif pairs
    combination_indices=get_permutations(df_motif.shape[0])
    df_mutation = pd.DataFrame(combination_indices, columns=['idx1', 'idx2'])
    df_mutation["region_idx"]=f"Region{region_idx}"
    df_mutation = df_mutation.merge(df_motif, left_on='idx1', right_index=True, suffixes=('', '1'))
    df_mutation = df_mutation.merge(df_motif, left_on='idx2', right_index=True, suffixes=('', '2'))
    df_mutation.rename(columns={
        'chromosome':'chromosome1',
        'start':'start1',
        'end':'end1',
        'start_rel':'start_rel1',
        'end_rel':'end_rel1',
        'protein':'protein1',
        "motif_length":"motif_length1"},inplace=True) 
    
    # remove overlapping motifs
    # if max(start1,start2)<=min(end1,end2), then the two motifs overlap, remove entire row
    idx_olap_row=df_mutation.apply(lambda row: max(row['start_rel1'],row['start_rel2'])<=min(row['end_rel1'],row['end_rel2']),axis=1)
    df_mutation=df_mutation[~idx_olap_row].reset_index(drop=True)
    if df_mutation.shape[0]==0:
        return
    
    # For each tf1 tf2 pair, get the null regions that is closest to tf1, record its location
    df_mutation.sort_values(by='start1',inplace=True)
    gr_tf1=pr.PyRanges(df_mutation.loc[:,['chromosome1','start1','end1']].rename(columns={"chromosome1":"Chromosome","start1":"Start","end1":"End"}))
    gr_tf1_null = gr_tf1.nearest(gr_null)
    df_mutation = pd.concat([df_mutation, gr_tf1_null.df.loc[:,["Start_b","End_b"]]], axis=1)
    df_mutation.rename(columns={"Start_b":"start_null","End_b":"end_null"},inplace=True)
    
    # for each null region longer than motif_length1, resize randomly to the same length as motif_length1
    df_mutation["max_shift_null"]=df_mutation.apply(lambda row: max(1,row["end_null"]-row["start_null"]-row["motif_length1"]-1),axis=1)
    # for each null region with max_shift_null>1, randomly shift the start position
    df_mutation.loc[df_mutation["max_shift_null"]>1,"start_null"]=df_mutation.loc[df_mutation["max_shift_null"]>1].apply(lambda row: np.random.randint(row["start_null"],row["end_null"]-row["motif_length1"]),axis=1)
    df_mutation.loc[df_mutation["max_shift_null"]>1,"end_null"]=df_mutation.loc[df_mutation["max_shift_null"]>1].apply(lambda row: row["start_null"]+row["motif_length1"],axis=1)   
    
    # calculate relative position of null region
    df_mutation["start_rel_null"]=df_mutation["start_null"]-region[1]
    df_mutation["end_rel_null"]=df_mutation["end_null"]-region[1]
    df_mutation.drop(columns=["max_shift_null","start_null","end_null"],inplace=True)
    
    # get sequences
    seq = seq_extractor.get_seq(region)
    df_mutation["seq_orig"]=seq
    df_mutation['seq_mut2'] = df_mutation.apply(lambda row: ablate_motifs(seq, row['start_rel2'], row['end_rel2']), axis=1)
    df_mutation['seq_mut1_null'] = df_mutation.apply(lambda row: ablate_motifs(seq, row['start_rel_null'], row['end_rel_null']), axis=1)
    df_mutation['seq_mut_both_null'] = df_mutation.apply(lambda row: ablate_motifs(seq, [row['start_rel_null'], row['start_rel2']],[row['end_rel_null'], row['end_rel2']]), axis=1)
    
    # get predictions
    df_pred_orig=get_predictions(df_mutation["seq_orig"],device=device,prefix="orig")
    df_pred_mut1=get_predictions(df_mutation["seq_mut1_null"],device=device,prefix="mut1_null")
    df_pred_mut2=get_predictions(df_mutation["seq_mut2"],device=device,prefix="mut2")
    df_pred_mut_both=get_predictions(df_mutation["seq_mut_both_null"],device=device,prefix="mut_both_null")
    
    # concatenate predictions with df_mutation by columns
    df_mutation=pd.concat([df_mutation,df_pred_orig,df_pred_mut1,df_pred_mut2,df_pred_mut_both],axis=1)

    # positive isa score indicate the original motif is contributing positively to RE activity.
    df_mutation=get_isas(df_mutation)
   
    df_mutation.drop(columns=['seq_orig',
                              'seq_mut2',
                              'seq_mut1_null',
                              'seq_mut_both_null'],inplace=True)
    # drop columns starting with pred
    df_mutation.drop(columns=[col for col in df_mutation.columns if col.startswith("pred")],inplace=True)
    # if out_path does not exist, create it using mode="w", include header
    if not os.path.exists(out_path):
        df_mutation.to_csv(out_path,mode="w",header=True,index=False)
        return
    df_mutation.to_csv(out_path,mode="a",header=False,index=False)
    return 




def full_analysis(file_prefix,device):
    # load data and tools
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    if "hepg2" in file_prefix:
        # chip file if needed: "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed"
        jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                        score_thresh=500,
                                        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv")
    elif "k562" in file_prefix:
        jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                        score_thresh=500,
                                        rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv")
    else:
        raise ValueError("file_prefix must contain either 'hepg2' or 'k562'")
    
    # read the regions
    df_regions = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_prefix}.bed",sep="\t",header=None)
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    for idx in range(df_regions.shape[0]):
        if idx%1000==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,[0,1,2]].tolist()
        write_mutation_res_for_one_region(seq_extractor,
                                          jaspar_annotator,
                                          region,
                                          f"df_mutate_pair_null_{file_prefix}.csv",
                                          idx,
                                          device)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_prefix",type=str)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    
    full_analysis(args.file_prefix,args.device)
    logger.info(f"Done with {args.file_prefix}")